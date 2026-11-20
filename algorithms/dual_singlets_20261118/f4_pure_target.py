"""Exact integral evaluation and a normalized target covector for the F4 fifth power."""
import argparse
import math
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat

import mixed_incoming as m
from f4_colour import basis
from f4_singlet_characters import singlets
from rank_four_s5 import decode_rows
from rank_four_sixth import G, IntegralFrame, pfaffian
from fortuity.hwcoh.evalg import table


@lru_cache(None)
def insertions(k,degree):
    left,right,sign,starts=table(k,degree,1)
    dest=np.repeat(np.arange(len(starts)),np.diff(np.r_[starts,len(left)]))
    return [(left[right==j],dest[right==j],sign[right==j]) for j in range(k)]


class SparseFrame(IntegralFrame):
    def mat_mul_letter(self, arr, degree, letter):
        result=np.zeros((1,math.comb(self.k,degree+1),27,27),np.int64)
        for (src,dst,sign),factors in zip(insertions(self.k,degree),self.factors[letter]):
            sub=arr[0,src]
            part=np.zeros_like(sub)
            for indices,values in factors:
                part+=np.take(sub,indices,axis=2)*values
            result[0,dst]+=part*sign[:,None,None]
        return result,degree+1


@lru_cache(None)
def matrices():
    mats,_=basis()
    result=[]
    for g in mats:
        a=g//np.gcd.reduce(np.abs(g).ravel())
        assert np.abs(a).sum(axis=1).max()<=2
        cols=[np.flatnonzero(c) for c in a.T]
        terms=[]
        for j in range(max(map(len,cols))):
            ids=np.asarray([x[j] if j<len(x) else 0 for x in cols])
            vals=np.asarray([a[x[j],i] if j<len(x) else 0 for i,x in enumerate(cols)])
            terms.append((ids,vals))
        result.append((a,terms))
    return result


def probe(task):
    words,alphabet,k,seed,target=task
    rng=np.random.default_rng(seed)
    fr=SparseFrame(27,1,k,[k],alphabet,seed=1)
    fr.factors={}
    # Coordinates in this frame are coefficients of the ordered top exterior product.
    palette=rng.choice(52,size=(9 if seed%2 else 18),replace=False)
    paired=rng.choice(palette,size=5,replace=False) if target is not None else None
    for letter in alphabet:
        select=rng.choice(palette,size=k)
        if paired is not None:
            select[0 if letter==G[0] else 1::2]=paired
        parts=[matrices()[i] for i in select]
        fr.L[letter]=(1,np.array([v for v,f in parts])[None])
        fr.factors[letter]=[f for v,f in parts]
    vals=fr.eval_words(words,dtype=np.int64)[:,0]
    if target is not None:
        sign,idx=target
        ab=np.einsum("iab,jba->ij",fr.L[G[0]][1][0],fr.L[G[1]][1][0])
        assert int(vals[idx])*sign==math.factorial(5)*pfaffian(ab-ab.T)
    return vals.tolist()


def evaluation(rows,k,first,n,workers,target):
    words=sorted({w for row in rows for w in row})
    alphabet=sorted({l for w in words for tr in w for l in tr})
    assert {m.nferm(w) for w in words}=={k} and {sum(map(len,w)) for w in words}=={k}
    bound=math.factorial(k)*27**(k//2)*2**k
    assert bound<2**63
    seeds=list(range(first,first+n))
    specification=(target[0],words.index(target[1])) if target else None
    with ProcessPoolExecutor(workers) as pool:
        vals=list(pool.map(probe,[(words,alphabet,k,s,specification) for s in seeds]))
    mat=fmpz_mat([[r.get(w,0) for w in words] for r in rows])*fmpz_mat(vals).transpose()
    return mat, bound,seeds


def covector(rows):
    exact=fmpq_mat(rows)
    b=exact.nrows()
    augmented=fmpq_mat([list(row)+[int(i==b-1)] for i,row in enumerate(exact.tolist())])
    reduced,rank=augmented.rref()
    sol=[0]*exact.ncols()
    for i in range(rank):
        j=next(j for j in range(reduced.ncols()) if reduced[i,j])
        assert j<exact.ncols()
        sol[j]=reduced[i,exact.ncols()]
    out=exact*fmpq_mat([[v] for v in sol])
    assert list(out.entries())==[0]*(b-1)+[1]
    return [[i,str(c)] for i,c in enumerate(sol) if c]


def run(out,workers):
    mu=(5,5,0,0,0)
    pd,sd=(singlets("F4",mu,n)[n] for n in (9,10))
    lo=m.prepare(mu,9)
    parents=decode_rows(lo,"sources")
    images=decode_rows(lo,"images")
    pure=decode_rows(m.prepare(mu,10),"sources")
    target=m.so_canon((G,)*5)
    sign,w=target
    pm,pbound,pseeds=evaluation(parents,9,325118,64,workers,None)
    all_rows=pure+images+[{w:sign}]
    full,bound,seeds=evaluation(all_rows,10,640088,96,workers,target)
    a=fmpz_mat(full.tolist()[:len(pure)])
    b=fmpz_mat(full.tolist()[len(pure):-1])
    c=fmpz_mat(full.tolist()[len(pure):])
    ranks=[pm.rank(),a.rank(),b.rank(),c.rank()]
    print("character",pd,sd,"evaluated ranks",ranks,flush=True)
    data=dict(weight=mu,lower_character=pd,terminal_character=sd,formal_lower=len(parents),
              formal_terminal=len(pure),lower_word_absolute_bound=pbound,terminal_word_absolute_bound=bound,
              lower_seeds=pseeds,terminal_seeds=seeds,ranks=ranks,
              lower_evaluations=[[int(x) for x in row] for row in pm.tolist()],
              terminal_evaluations=[[int(x) for x in row] for row in full.tolist()],
              fifth_power_pfaffian_checks=len(seeds))
    if ranks==[pd,sd,sd-1,sd]:
        data["normalized_target_covector"]=covector(c)
    m.write(out,data)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    assert __debug__
    run(**vars(p.parse_args()))
