"""Replay one fixed level-34 first correction on fresh and on full coupled profile frames.

Each coefficient multiplies its original projected trace recipe *before* taking Q0.
Residual orientation is (-Q1 O0) - Q0 O1. No solve is performed by this program.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
import projected_q0 as maps

HERE=Path(__file__).resolve().parent
PRIME=1000003
PLAN_SHA='1b7ac3dec77fcb50473d01176a8f309773e8411f0fd4e969f58299a026df6e8c'
MINOR_SHA='6db12024f5a86596a8c161628e0661097dc3b11b9cbee49019e9b9ce9542f199'
SOLUTION_SHA='cfeeb05b32766abf796d7fad1b7f08cf136a3b60746fda598d647c31c0305ca2'
RHS_SHA='9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14'
TARGETS=(13,14,15,16)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def vector_sha(v):return hashlib.sha256(np.asarray(v,dtype='<u4').tobytes()).hexdigest()

@lru_cache(maxsize=512)
def coeff(c):
    r=Fraction(c);return (r.numerator*pow(r.denominator,-1,PRIME))%PRIME

def geometry(plan,src,t):
    from g2_full_classical_screen_v1 import contents, sector
    charge=sector(tuple(plan['mu']),plan['source_length'])
    alphabet=sorted({l for n in (plan['source_length'],plan['target_length']) for c in contents(tuple(plan['mu']),n) for l in c})
    assert [i for i in range(len(plan['targets'])) if any(i in plan['sources'][s]['targets'] for s in src)] == list(TARGETS)
    return charge,alphabet

def make_frame(t,alphabet,points,seed):
    from fortuity.full.evalg import Frame
    g=t['fermions']+2
    f=Frame(3,points,g,[t['fermions']],alphabet,seed=seed,p=PRIME)
    # The pair of bosonic matrices fails to commute at each point.
    matrices=[v[:,0] for d,v in f.L.values() if d==0];a,b=matrices[:2]
    assert np.all(np.any((a@b-b@a)%PRIME,axis=(1,2)))
    return f

def calculate(args):
    sys.path.insert(0,str(args.runtime.resolve()))
    import su_trace_adapter_v1 as algebra
    from native_su_parallel_v1 import evaluate, evaluate_profile
    from native_su_projection_v1 import project
    from sp_full_profiles import profile
    raw=args.plan;assert sha(raw)==PLAN_SHA and sha(args.minor)==MINOR_SHA and sha(args.solution)==SOLUTION_SHA
    plan=json.loads(raw.read_bytes())['incoming'];m=json.loads(gzip.decompress(args.minor.read_bytes()));s=json.loads(gzip.decompress(args.solution.read_bytes()))
    src=m['map']['sources'];assert src==[19,20,21,22] and m['map']['targets']==[16] and m['determinant']==258962
    assert s['base_minor_sha256']==MINOR_SHA and s['profile']==16 and s['prime']==m['map']['prime']==PRIME
    assert len(s['base_coefficients'])==len(m['rows'])==7211 and not s['extra_rows'] and not s['extra_coefficients']
    t=plan['targets'][args.target]; assert t['index']==args.target
    charge,alphabet=geometry(plan,src,t);pattern=maps.pattern(t)
    seed=2026111801+7919*t['index'] if args.mode=='production' else 88304773+229*t['index']
    points=math.ceil((t['dimension']+8)/math.comb(t['fermions']+2,t['fermions'])) if args.mode=='production' else args.points
    f=make_frame(t,alphabet,points,seed);w=f.ncols
    # A direct short prefix on one point is evaluated through the independent Python
    # exterior algebra. Use the identical letter frame, not a different RNG prefix.
    from copy import copy
    first=copy(f);first.P=1;first.ncols=first.cpp;first.L={l:(d,v[:1]) for l,(d,v) in f.L.items()}
    h=hashlib.sha256();total=terms=0;buf={};rhs=np.zeros(w,dtype=np.int64); controls=0
    with gzip.open(args.rhs,'rb') as inp:
        for line in inp:
            total+=1;h.update(line);c,raw_word=json.loads(line);word=maps.decode(raw_word)
            if profile(tuple(l for tr in word for l in tr))!=pattern:continue
            assert word not in buf
            buf[word]=coeff(c);terms+=1
            if controls<3:
                vv=first.eval_words([word],dtype=np.int64)[0]*coeff(c)%PRIME
                assert np.array_equal(vv,evaluate([{word:coeff(c)}],first,1)[0]);controls+=1
            if len(buf)==256:
                rhs=(rhs+evaluate([buf],f,args.workers)[0])%PRIME;buf={}
    if buf:rhs=(rhs+evaluate([buf],f,args.workers)[0])%PRIME
    assert total==517468 and terms==plan['first_rhs_trace_terms'][str(t['index'])] and h.hexdigest()==RHS_SHA and controls==3
    lhs=np.zeros(w,dtype=np.int64);used=0;checks=0;contrib=Counter()
    rows=list(zip(m['rows'],s['base_coefficients']))
    for start in range(0,len(rows),32):
        expressions={}
        for (rel,raw_word),c in rows[start:start+32]:
            absolute=src[rel];sr=plan['sources'][absolute]
            if t['index'] not in sr['targets']:continue
            assert c and 0<c<PRIME;word=maps.decode(raw_word)
            exp=project(algebra.parity_project({word:1},1),charge,PRIME)
            assert exp;used+=1;contrib[absolute]+=1
            if not checks:
                direct=algebra.differential(exp);sub={v:n for v,n in direct.items() if profile(tuple(l for tr in v for l in tr))==pattern}
                if sub:
                    # Both sides on an identical one-point frame; the direct path uses
                    # the Python graded word differential before trace evaluation.
                    expected=evaluate([sub],first,1)[0];actual=evaluate_profile([exp],first,maps.pattern(sr),pattern,1)[0]
                    assert np.array_equal(expected,actual)
                    if actual.any():checks+=1
            out=expressions.setdefault(absolute,Counter())
            for v,n in exp.items():out[v]+=int(c)*int(n)
        for absolute,out in expressions.items():
            combined={v:n%PRIME for v,n in out.items() if n%PRIME}
            if combined:lhs=(lhs+evaluate_profile([combined],f,maps.pattern(plan['sources'][absolute]),pattern,args.workers)[0])%PRIME
    assert checks==1 and used==sum(contrib.values())
    residual=(rhs-lhs)%PRIME;nonzero=np.flatnonzero(residual); sample=None if not len(nonzero) else int(nonzero[0]);stem=f'{args.mode}_{args.target}'
    args.out.mkdir(exist_ok=True,parents=True)
    saved=np.vstack((rhs,lhs,residual)).astype('<u4');np.save(args.out/(stem+'.npy'),saved,allow_pickle=False)
    if args.target==16 and args.mode=='production' and args.production_rhs is not None:
        assert np.array_equal(rhs,np.load(args.production_rhs,allow_pickle=False))
    record=dict(target=args.target,mode=args.mode,prime=PRIME,seed=seed,noncommuting_points=points,exterior_coordinates_per_point=f.cpp,
        width=w,target_dimension=t['dimension'],rhs_terms=terms,source_recipes_used=used,recipes_by_source=dict(sorted(contrib.items())),
        solution_sha256=SOLUTION_SHA,source_minor_sha256=MINOR_SHA,plan_sha256=PLAN_SHA,formal_rhs_stream_sha256=RHS_SHA,
        rhs_sha256=vector_sha(rhs),lhs_sha256=vector_sha(lhs),residual_sha256=vector_sha(residual),rhs_nonzero=int(np.count_nonzero(rhs)),
        lhs_nonzero=int(np.count_nonzero(lhs)),residual_nonzero=len(nonzero),first_nonzero=None if sample is None else dict(index=sample,rhs=int(rhs[sample]),lhs=int(lhs[sample]),residual=int(residual[sample])),
        direct_controls=dict(rhs_words=controls,projected_source=checks),vectors_sha256=sha(args.out/(stem+'.npy')))
    (args.out/(stem+'.json')).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--target',type=int,choices=TARGETS,required=True);p.add_argument('--mode',choices=('fresh','production'),default='fresh');p.add_argument('--points',type=int,default=12);p.add_argument('--workers',type=int,default=4)
    for name,path in [('plan','classical_order1.json'),('rhs','negative_Q1.jsonl.gz'),('minor','other_maps/task_3/minor.json.gz'),('solution','first_profile16/solution.json.gz')]:p.add_argument('--'+name,type=Path,default=HERE/path)
    p.add_argument('--production-rhs',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert a.points>0 and a.workers>0;calculate(a)
