"""Evaluate the first-equation residual on integer traceless matrices killing selected parent images.

Every substitution is made in the matrix entries before taking traces.  Making
both letters in a split diagonal kills their graded commutator, including in an
exterior coefficient algebra. All words in that parent content then vanish on
this target under the projected Q0, whatever their rational coefficients.
"""
import argparse, hashlib, json
from collections import Counter
from fractions import Fraction
from pathlib import Path
from random import Random
import finite_trace_identity as ft
import formal_trace_blocks as f
import replay_formal_blocks as s
from fortuity.full.words import Q_multi

TARGETS=('0000000000010204060c16','0000000000010204060e14')

def split_pairs(t):
    t=bytes.fromhex(t); parents={}
    for i in range(len(t)):
        for j in range(i+1,len(t)):
            a,b=t[i],t[j]; k=f.MERGE.get(bytes((a,b)))
            if k is None:continue
            g=bytes(sorted([v for h,v in enumerate(t) if h not in (i,j)]+[k]));terms={tuple(pair):v for v,pair in f.SPLITS[k] if sorted(pair)==[a,b]}; assert a!=b and set(terms)=={(a,b),(b,a)} and terms[b,a]==-terms[a,b]*(-1)**(ft.ODD[a]*ft.ODD[b]); info=(k,a,b,terms[a,b]);assert parents.setdefault(g.hex(),info)==info
    return parents

def vector(shared):
    z=s.read(shared);D=int(z['denominator']);v=[(s.unpack_hex(w),Fraction(int(c),D)) for w,c in z['words']];assert D>0 and len(v)==len(set(w for w,c in v)) and all(c for w,c in v);vhash=hashlib.sha256(json.dumps([[s.hexes(w),str(c)] for w,c in v],separators=(',',':')).encode()).hexdigest();return v,vhash

def residuals(rhs,shared,targets):
    dom=set(map(bytes.fromhex,targets));rows={g:{} for g in dom};n=0;h=hashlib.sha256();v,hv=vector(shared)
    def content(w):return bytes(sorted(b''.join(w)))
    for w,c in s.stream_rows(rhs,h):
        n+=1;g=content(w)
        if g in rows:assert w not in rows[g];rows[g][w]=-Fraction(c)
    assert n==517468 and h.hexdigest()==s.INPUT
    needed={bytes.fromhex(k) for t in targets for k in split_pairs(t)}
    for w,c in v:
        if content(w) not in needed:continue
        for u,a in Q_multi(s.expanded(w)).items():
            z=s.packed(u);g=content(z)
            if g in rows:rows[g][z]=rows[g].get(z,0)+a*c
    return hv,{g.hex():{w:c for w,c in r.items() if c} for g,r in rows.items()}

def integer_test(t,row,parents,allowed,seed):
    active={i for g in allowed for i in parents[g][1:3]};ls=sorted(set(bytes.fromhex(t)));n=sum(ft.ODD[i] for i in bytes.fromhex(t));assert n<=3,'factorial integer control is intended for small exterior degrees';rnd=Random(seed);M={}
    for l in ls:
        for h in (range(n) if ft.ODD[l] else [-1]):
            a=ft.traceless(rnd,3); M[l,h]=tuple(tuple(x if i==j or l not in active else 0 for j,x in enumerate(z)) for i,z in enumerate(a))
    # Check each homogeneous coefficient matrix commutes in the ordinary matrix algebra.
    # Moving the exterior coefficients past one another supplies the graded sign.
    for g in allowed:
        k,a,b,co=parents[g]
        for (l,h),A in M.items():
            if l!=a:continue
            for (ll,hh),B in M.items():
                if ll==b:assert ft.mm(A,B)==ft.mm(B,A)
    value=sum((c*ft.exterior_value(w,M,n) for w,c in row.items()),Fraction());return dict(seed=seed,exterior_degree=n,diagonal_letters=bytes(sorted(active)).hex(),matrices={f'{l:02x}:{h}':v for (l,h),v in M.items()},rational_value=str(value))

def main(rhs,shared,frontier):
    z=s.read(frontier);v,rr=residuals(rhs,shared,TARGETS);assert v==z['vector_sha'] and s.INPUT==z['rhs'];out={}
    for t in TARGETS:
        p=split_pairs(t);entry=z['target_parents'][t];assert sorted(p)==entry['all'];safe=entry['confined_to_failures'];assert len(safe)==1; controls={}
        for label,allowed in [('confined',safe),('confined_plus_one_success',safe+entry['touch_one_success']),('all',sorted(p))]:
            controls[label]=dict(parents=allowed,assignments=[integer_test(t,rr[t],p,allowed,seed) for seed in (817234,419071)])
        out[t]=dict(formal_residual_words=len(rr[t]),splits_by_parent={k:[f'{i:02x}' for i in r[:3]]+[r[3]] for k,r in sorted(p.items())},tests=controls)
    return dict(rhs=s.INPUT,vector_sha=v,targets=out,sign='Q0 O1 - (-Q1 O0)',split_order=['merged_letter','left','right','coefficient_of_left_right'],scope='Nonzero integer-specialization value excludes the tested parent subspace on this target. A zero value is inconclusive.')

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--shared',type=Path,default=s.HERE/'joint_70/vector_shared.json.gz');a.add_argument('--frontier',type=Path,default=s.HERE/'joint_70/parent_frontier.json');a.add_argument('--out',type=Path,required=True);x=a.parse_args();assert __debug__;r=main(x.rhs,x.shared,x.frontier);s.write(x.out,r);print({t:{k:[a['rational_value'] for a in z['assignments']] for k,z in d['tests'].items()} for t,d in r['targets'].items()})
