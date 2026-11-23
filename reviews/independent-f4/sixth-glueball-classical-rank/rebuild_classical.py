"""Rebuild selected Albert contractions with integer matrices and sparse tensor actions."""
from functools import cache
from itertools import combinations,permutations,product
from pathlib import Path
import gzip
import json

import numpy as np
from scipy.sparse import csr_matrix

from colour_tensors import act,invariant_ends
from exceptional_colour import projectors
from double_trace_flag import Evaluation
from triple_trace_flag import Triple
from single_trace_flag import native26

HERE=Path(__file__).resolve().parent
FOUR_END=((0,1,2),(0,2,3),(1,2,0),(2,0,1),(3,0,2))
FOUR_MID=((0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,0,1),(0,1,0))
SIX_END=((0,1,2),(0,2,3),(2,0,1),(1,0,3))
SIX_MID=((0,1,2),(0,1,3),(1,2,3),(0,0,1))
SIX_PAIRS=((0,1),(0,2),(1,3))


def sign(seq):
    return (-1)**sum(a>b for k,a in enumerate(seq) for b in seq[k+1:])


def flags(row):
    stage,i=row
    if stage=="four":
        return None,tuple(product(FOUR_END,FOUR_MID,FOUR_END,FOUR_MID))[i]
    assert stage=="nested"
    return SIX_PAIRS[i//256],tuple(product(SIX_END,SIX_MID,SIX_END,SIX_MID))[i%256]


def trace_frame(rows):
    catalogue=json.loads((HERE/"triple_selected.json").read_bytes())["rows"]
    chosen=[catalogue[i] for s,i in rows if s=="reference"]
    ordinary=[v for v in chosen if not v[0].startswith("T")]
    triples=[v for v in chosen if v[0].startswith("T")]
    first,second=Evaluation(ordinary),Triple(triples)
    at={tuple(x):i for i,x in enumerate(ordinary+triples)}
    ix=[at[tuple(x)] for x in chosen]
    return lambda n,p:np.r_[first.point(n,p),second.point(n,p)][ix]


def graph_point(rep,ends,plan,n,p):
    low,high,exchange=ends
    coeff=np.random.Generator(np.random.PCG64(144336+n)).integers(-3,4,size=(4,52))
    x=(coeff@rep.reshape(52,-1)).reshape(4,26,26)%p
    marks={(i,j):(x[i]@x[j]-x[j]@x[i])%p for i,j in combinations(range(4),2)}
    def oriented(mat,side):return csr_matrix(mat.T if side else mat)
    xs=[[oriented(m,s) for m in x] for s in range(2)]
    ms=[{k:oriented(m,s) for k,m in marks.items()} for s in range(2)]
    @cache
    def endpoint(side,arms):
        total=np.zeros_like(low)
        for i,j in combinations(range(4),2):
            k,l=(r for r in range(4) if r not in (i,j))
            u=act(ms[side][i,j],low if side else high,arms[0],p)
            a=act(xs[side][l],act(xs[side][k],u,arms[1],p),arms[2],p)
            b=act(xs[side][k],act(xs[side][l],u,arms[1],p),arms[2],p)
            total+=sign((i,j,k,l))*(a-b)
        return total%p
    @cache
    def half(side,end,arms):
        total=np.zeros_like(low)
        for per in permutations(range(3)):
            u=endpoint(side,end)
            for k,axis in zip(per,arms):u=act(xs[side][k],u,axis,p)
            total+=sign(per)*u
        return (total%p).reshape(-1)
    @cache
    def right(end,mid,axes):
        val=half(False,end,mid)
        if axes is None:return val
        return act(exchange,val.reshape((26,)*4),axes,p).ravel()
    assert 26**4*(p-1)**2<2**63
    result=np.array([int(half(True,l,lm)@right(r,rm,axes))%p
                     for axes,(l,lm,r,rm) in plan])
    endpoint.cache_clear();half.cache_clear();right.cache_clear()
    return result


def regenerate(p,info,positions):
    rows=info["rows"]
    ref=trace_frame(rows)
    plan=[flags(r) for r in rows if r[0]!="reference"]
    rep=native26().astype(np.int64)
    ends=invariant_ends(rep,*projectors(HERE,rep,p)[12],p)
    for pos in positions:
        assert 0<=pos<703
        n=info["points"][pos]
        yield pos,np.r_[ref(n,p),graph_point(rep,ends,plan,n,p)]
