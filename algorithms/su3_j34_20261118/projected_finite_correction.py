"""Exact projection of the first correction to specified complete letter contents.

Reduce products of four or more traces by the graded three-colour identity.
Retain the identity rows used as well as the original coefficient space. A
nonzero remainder in this presentation is not a physical separating witness.
"""
import argparse, gzip, hashlib, json
from fractions import Fraction
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path
import replay_formal_blocks as src
from fortuity.full.words import Q_multi
import finite_trace_identity as ft
HERE=Path(__file__).resolve().parent
FIRST=('0000000000000505050711','0000000000000505050814','000000000000050505090f','0000000000000505050c10')

def content(w):return bytes(sorted(b''.join(w)))

def clean(x):return {k:v for k,v in x.items() if v}

def axpy(dst,row,c):
    for k,v in row.items():
        n=dst.get(k,0)+c*v
        if n:dst[k]=n
        else:dst.pop(k,None)

@lru_cache(100000)
def replacement(w):
    if len(w)<4:return {w:1}
    rel=ft.identity(w[:4]);assert rel.pop(w[:4])==1
    r={}
    for u,v in rel.items():
        s,t=ft.multi(u+w[4:]);
        if s:axpy(r,replacement(t),-s*v)
    assert all(len(t)<=3 for t in r);return r

def reduce(row):
    d={}
    for w,c in row.items():axpy(d,replacement(w),c)
    return d

def rhs_rows(rhs,contents):
    d={};h=hashlib.sha256();n=0
    for w,c in src.stream_rows(rhs,h):
        n+=1
        if content(w) in contents:
            assert w not in d;d[w]=Fraction(c)
    assert n==517468 and h.hexdigest()==src.INPUT;return d

def project(w,contents):
    return {t:c for v,c in Q_multi(src.expanded(w)).items() if c and content(t:=src.packed(v)) in contents}

def rref(images,rhs):
    """Exact sparse column elimination, with solutions as original column coordinates."""
    try: from flint import fmpq
    except ImportError: fmpq=Fraction
    def q(x):return fmpq(x.numerator,x.denominator) if isinstance(x,Fraction) else fmpq(x)
    # Use frequent coordinates later to reduce fill; ordering is independent of RHS values.
    freq=Counter(k for r in images for k in r);coords=sorted(freq.keys() | rhs.keys(),key=lambda k:(freq[k],k));idx={k:i for i,k in enumerate(coords)}
    piv={};unit=fmpq(1);coeffs=[]
    def clear(y,cs):
        while y:
            i=min(y)
            if i not in piv:break
            row,position,prior=piv[i];v=y[i];axpy(y,row,-v);cs[position]=cs.get(position,0)+v
        return y,cs
    for j,img in enumerate(images):
        y={idx[k]:q(v) for k,v in img.items()};prior={};clear(y,prior)
        if y:
            i=min(y);d=y[i];piv[i]=({k:v/d for k,v in y.items()},len(coeffs),prior);coeffs.append((j,d,prior))
    y={idx[k]:q(v) for k,v in rhs.items()}; sol={};clear(y,sol)
    result={};res={coords[i]:str(v) for i,v in y.items()}
    for i in reversed(range(len(coeffs))):
        val=sol.get(i,0)
        if not val:continue
        j,d,prior=coeffs[i];v=val/d;result[j]=Fraction(str(v));axpy(sol,prior,-v)
    return dict(rank=len(piv),coordinates=len(coords),augmented_rank=len(piv)+bool(y),residual_terms=len(y)),result,res

def correction(rhs,contents):
    d=rhs_rows(rhs,contents); parents=sorted(set().union(*(src.parents(src.expanded(w)) for w in d))); print('parents',len(parents),flush=True)
    images=[project(w,contents) for w in parents]; raw,_,_=rref(images,d); red,coef,rem=rref([reduce(i) for i in images],reduce(d));print('formal',raw,'identity-reduced',red,flush=True)
    diff=d.copy()
    for i,v in coef.items():axpy(diff,images[i],-v)
    assert {k:Fraction(v) for k,v in rem.items()}==reduce(diff)
    rows=[[src.hexes(parents[i]),str(c)] for i,c in sorted(coef.items())]
    return dict(rhs_sha256=src.INPUT,contents=[i.hex() for i in sorted(contents)],projected_terms=len(d),parents=len(parents),formal=raw,
                three_trace=red,solution_or_partial=rows,formal_residual_before_relations=len(diff),after_relations=len(rem),
                other_contents_enforced=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rhs',type=Path,default=HERE/'negative_Q1.jsonl.gz');p.add_argument('--contents',nargs='+',default=FIRST);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert __debug__
    r=correction(a.rhs,set(map(bytes.fromhex,a.contents)));src.write(a.out,r)
