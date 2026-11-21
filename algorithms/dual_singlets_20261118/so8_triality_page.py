"""Classical highest-weight bounds by triality for a fixed SO(8) derivative grade."""
import argparse
import math
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from fortuity.full.letters import is_ferm
from so8_triality import representations
from trace_descent import dimensions, highest_weights, multiply_rows


def probe(rows, sizes, seed, spinors):
    parts=[[],[],[]]
    for degree,dimension in sizes.items():
        selected=[{w:c for w,c in row.items() if m.nferm(w)==degree} for row in rows]
        words=sorted({w for row in selected for w in row})
        if not words:
            continue
        alphabet=sorted({l for w in words for tr in w for l in tr})
        # k=d+1 gives d+1 exterior coordinates per point.
        k=degree+1
        assert 8*(m.P-1)**2 < 2**63 and math.comb(k,k//2)*(m.P-1)**2<2**63
        points=math.ceil((dimension+2*k+8)/k)
        for start in range(0,points,5):
            count=min(5,points-start)
            rng=np.random.default_rng(seed+6007*start+8191*degree)
            frame=m.Frame(8,count,k,[degree],[],p=m.P)
            vector={}
            up=np.triu_indices(8,1)
            for l in alphabet:
                odd=int(is_ferm(l))
                a=np.triu(rng.integers(0,m.P,(count,k if odd else 1,8,8),dtype=np.int64),1)
                vector[l]=odd,(a-a.swapaxes(-1,-2))%m.P
            for i,spin in enumerate((None,*spinors)):
                frame.L=vector if spin is None else {
                    l:(odd,np.tensordot(a[...,up[0],up[1]],spin,axes=1)*pow(2,-1,m.P)%m.P)
                    for l,(odd,a) in vector.items()}
                result=multiply_rows(selected,words,frame.eval_words(words,dtype=np.int64))
                parts[i].append(np.asarray(result%m.P,dtype=np.int64))
    if not parts[0]:
        return [np.zeros((len(rows),0),dtype=np.int64) for _ in range(3)]
    vector,left,right=[np.column_stack(x) for x in parts]
    return [vector,(vector+left+right)%m.P,(2*vector-left-right)%m.P]


def all_rows(weight,length):
    return [row for block in highest_weights(weight,length).values() for row in block]


def rank(a):
    return nmod_mat(a.tolist(),m.P).rank() if a.size else 0


def run(grade,deficit,seed,*,taylor_weight=None):
    w=(grade,grade)+(4-grade,)*3 if taylor_weight is None else tuple(taylor_weight)
    if len(w) != 5 or min(w) < 0 or w[0] < w[1] or list(w[2:]) != sorted(w[2:],reverse=True):
        raise ValueError("expected dominant spin and flavour degrees")
    y=sum(w)-deficit
    source=all_rows(w,y-1)
    target=all_rows(w,y)
    following=dimensions(w,y+1,8,"SO")
    target_sizes=dimensions(w,y,8,"SO")
    spinors=representations()
    incoming=[m.q(row) for row in source]
    outgoing=[m.q(row) for row in target]
    assert all(not m.q(row) for row in incoming+outgoing)
    measured=probe(target+incoming,target_sizes,seed,spinors)
    after=probe(outgoing,following,seed+14231,spinors)
    res={}
    for name,mat,image in zip(("even","invariant","standard_even"),measured,after,strict=True):
        cochains=rank(mat[:len(target)])
        into=rank(mat[len(target):])
        out=rank(image)
        assert cochains>=into+out
        res[name]=dict(cochains=cochains,incoming_minor_rank=into,outgoing_minor_rank=out,
                       classical_cohomology_upper=cochains-into-out)
    assert res["even"]["cochains"]==sum(target_sizes.values())
    assert res["even"]["cochains"]==res["invariant"]["cochains"]+res["standard_even"]["cochains"]
    return dict(weight=w,length=y,target_deficit=deficit,grade=grade,prime=m.P,seed=seed,
                target_character=target_sizes,following_character=following,
                formal_preceding=len(source),formal_target=len(target),ranks=res)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grade",required=True,type=int,choices=range(5))
    p.add_argument("--deficit",type=int,required=True,choices=range(5))
    p.add_argument("--seed",type=int,default=241093)
    p.add_argument("--taylor-weight",type=int,nargs=5,help="evaluate the specified Taylor weight at the same deficit")
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert __debug__
    result=run(a.grade,a.deficit,a.seed,taylor_weight=a.taylor_weight)
    m.write(a.out,result)
    print(result["ranks"],flush=True)
