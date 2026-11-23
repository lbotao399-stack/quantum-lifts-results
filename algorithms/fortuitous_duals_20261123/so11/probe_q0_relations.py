#!/usr/bin/env python3
"""Exact coefficients of Q0 images of genuine SO11 Gram identities.

Every source is a degree14 Gram determinant with twelve marked endpoints on
each side and one matter jet psi_BC. Q0 replaces that edge by BC-CB. Scalar
projection and six odd returns commute with Q0 and finite-colour restriction.
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
import argparse,hashlib,json,random,subprocess,time
import independent_trace_ops as tr

HERE=Path(__file__).resolve().parent

def source(index):
    rng=random.Random(971235+104729*index)
    species=list('AAAAABBBBCCCC')+['x']
    # Random determinant endpoint ownership, retaining two ordinary edges.
    endpoints=list(range(28));rng.shuffle(endpoints)
    return dict(index=index,species=species,endpoints=endpoints,
                matter_jet='x has Q0(x)=BC-CB',determinant_size=12)

def q_inputs(spec,mincycles):
    u=[v if v!=27 else 29 for v in spec['endpoints']]
    u+= [27,28]
    words=[spec['species'][:-1]+list('BC'),spec['species'][:-1]+list('CB')]
    for word,sg in zip(words,[1,-1]):
        s=f'{mincycles} 0\n'+' '.join(str('ABC'.index(c)) for c in word)+'\n'
        s+=' '.join(['0']*15)+'\n'+' '.join(map(str,u))+'\n'
        yield sg,s

def run(spec,mincycles=6):
    blocks={tuple(sorted(map(len,q))):defaultdict(Fraction) for q in tr.QUERIES}
    counters=[];hashes=[]
    for sg,inp in q_inputs(spec,mincycles):
        raw=subprocess.check_output([str(HERE/'probe_trials/bmngram')],input=inp,text=True)
        hashes.append(hashlib.sha256(raw.encode()).hexdigest())
        for line in raw.splitlines():
            traces,coeff=line.split()
            if traces=='!':counters.append(int(coeff));continue
            words=tuple(traces.rstrip('.').split('.'))
            shape=tuple(sorted(map(len,words)))
            if shape not in blocks:continue
            s,w=tr.canonical(words)
            if s:blocks[shape][w]+=sg*s*int(coeff)
    answers=[];sizes=[]
    for query in tr.QUERIES:
        shape=tuple(sorted(map(len,query)))
        poly={w:c for w,c in blocks[shape].items() if c}
        sizes.append(len(poly));poly=tr.project(poly)
        for scalar in 'ABC':
            for ghost in '-+':poly=tr.action(poly,scalar,ghost)
        sign,w=tr.canonical(query)
        answers.append(str(sign*poly.get(w,0)))
    return dict(source=spec,query_values=answers,canonical_raw_terms=sizes,
                matching_counters=counters,raw_sha256=hashes)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--count',type=int,default=24);ap.add_argument('--start',type=int,default=0)
    args=ap.parse_args();results=[];t=time.monotonic()
    for i in range(args.start,args.start+args.count):
        s=time.monotonic();r=run(source(i));r['seconds']=time.monotonic()-s;results.append(r)
        (HERE/'probe_q0_relations.json').write_text(json.dumps({'rows':results,'seconds':time.monotonic()-t},indent=2)+'\n')
        print(i,r['query_values'],r['matching_counters'],r['seconds'],flush=True)
if __name__=='__main__':main()
