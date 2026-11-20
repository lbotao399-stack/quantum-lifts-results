"""Terminal spin/flavour singlets and rank-four symplectic Pfaffian relations at level 30."""
import argparse
from collections import Counter
from functools import cache
from itertools import permutations
from math import gcd
from pathlib import Path
import random

from flint import fmpz_mat

import mixed_incoming as m
from orthogonal_relations import lower
from symplectic_relations import matchings
from terminal_commutators import parity, sign, words
from replay_orthogonal import exchanges


def bounds(t):
    weight=(t,t)+(5-t,)*3
    formal=cartan=0
    for spin,ss in ((0,1),(1,-1)):
        for perm in permutations(range(3)):
            mu=(t+spin,t-spin)+tuple(5-t+perm[i]-i for i in range(3))
            row=words(mu)
            c=ss*sign(perm)
            formal+=c*len(row)
            cartan+=c*sum(len(w)<=4 for w in row)
    return dict(t=t,weight=weight,length=sum(weight),formal=formal,cartan=cartan)


def quotient(weight,dim):
    basis=words(weight)
    ids={w:i for i,w in enumerate(basis)}
    rows=[]
    for i,j in ((0,1),(2,3),(3,4)):
        shifted=list(weight);shifted[i]+=1;shifted[j]-=1
        for w in words(tuple(shifted)):
            coeff=lower(w,i,j)
            row=[0]*len(basis)
            for v,c in coeff.items():
                if c:
                    row[ids[v]]=c
            rows.append(row)
    matrix=fmpz_mat(rows)
    kernel,k=matrix.nullspace()
    assert k==dim and matrix.rank()+k==len(basis)
    cols=[]
    for j in range(k):
        v=[int(kernel[i,j]) for i in range(len(basis))]
        den=gcd(*v)
        cols.append([x//den for x in v])
    cov=fmpz_mat(cols).transpose()
    assert (matrix*cov).is_zero()
    return basis,cov


@cache
def pf_terms():
    return tuple(matchings(tuple(range(10))))


def relation(species,slots):
    n=len(species)
    assert sorted(slots)==list(range(2*n)) and len(set(i//2 for i in slots[:10]))==10
    base_mate=[0]*(2*n)
    base_sign=[1]*(2*n)
    for a,b in zip(slots[10::2],slots[11::2]):
        base_mate[a],base_mate[b]=b,a
        base_sign[b]=-1
    ans=Counter()
    for pairs,c in pf_terms():
        mate=base_mate.copy()
        edge=base_sign.copy()
        for i,j in pairs:
            a,b=slots[i],slots[j]
            mate[a],mate[b]=b,a
            edge[b]=-1
        visited=0
        terms=[]
        for first in range(n):
            if visited&(1<<first):
                continue
            start=pos=2*first
            counts=[0]*5
            odd=[]
            while True:
                letter=pos//2
                visited|=1<<letter
                species_id=species[letter]
                counts[species_id]+=1
                if species_id<2:
                    odd.append(letter)
                c*=edge[pos^1]
                pos=mate[pos^1]
                if pos==start:
                    break
            if sum(counts)%2 or max(counts[:2])>1:
                break
            terms.append((tuple(counts),sorted(odd,key=species.__getitem__)))
        else:
            terms.sort(key=lambda block:(sum(block[0]),block[0]))
            word=tuple(block[0] for block in terms)
            if any(a==b and parity(a) for a,b in zip(word,word[1:])):
                continue
            ans[word]+=c*(-1)**exchanges([i for _,labels in terms for i in labels])
    return {w:c for w,c in ans.items() if c}


def make_case(weight,rng,trial):
    species=[i for i,n in enumerate(weight) for _ in range(n)]
    rng.shuffle(species)
    order=list(range(len(species)))
    rng.shuffle(order)
    free=[2*i+rng.randrange(2) for i in order[:10]]
    rng.shuffle(free)
    fixed=[i^1 for i in free]
    for i in order[10:]:
        fixed.extend((2*i,2*i+1))
    if trial%3:
        rng.shuffle(fixed)
    return species,free+fixed


def run(out,attempts,seed):
    screen=[bounds(t) for t in range(6)]
    case=screen[3]
    weight=case["weight"]
    basis,cov=quotient(weight,case["formal"])
    ids={w:i for i,w in enumerate(basis)}
    rng=random.Random(seed)
    selected=[]
    projected=[]
    for t in range(attempts):
        species,slots=make_case(weight,rng,t)
        coeff=relation(species,slots)
        if not coeff:
            continue
        raw=[0]*len(basis)
        for w,c in coeff.items():
            raw[ids[w]]=c
        row=list(map(int,(fmpz_mat([raw])*cov).tolist()[0]))
        if any(row) and fmpz_mat(projected+[row]).rank()>len(projected):
            projected.append(row)
            selected.append(dict(trial=t,species=species,slots=slots,
                                 raw=[[i,c] for i,c in enumerate(raw) if c],singlet=row))
            print("symplectic relation",len(selected),"trial",t,flush=True)
            if len(projected)==case["formal"]-case["cartan"]:
                break
    assert len(projected)<=case["formal"]-case["cartan"]
    minor={}
    if projected:
        from itertools import combinations
        columns,det=next((cols,int(sq.det())) for cols in combinations(range(cov.ncols()),len(projected))
                         for sq in [fmpz_mat([[r[i] for i in cols] for r in projected])] if sq.det())
        minor=dict(columns=columns,determinant=det)
    m.write(out,dict(level=30,terminal=screen,group="USp(8)",weight=weight,weight_dimension=len(basis),
                     pfaffian_terms=len(pf_terms()),attempts=t+1,seed=seed,relations=selected,minor=minor,
                     classical_lower=case["cartan"],classical_upper=case["formal"]-len(projected)))
    print("mixed terminal cohomology bounds",case["cartan"],case["formal"]-len(projected),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--attempts",type=int,default=1600)
    p.add_argument("--seed",type=int,default=1627924)
    assert __debug__
    run(**vars(p.parse_args()))
