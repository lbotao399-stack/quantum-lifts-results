"""Reflection-odd terminal Pfaffians and both SO(8) parities at level 30."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from math import prod
from pathlib import Path

import numpy as np

import mixed_incoming as m
from mixed33_exact import read
from rank_four_sixth import pfaffian
from rank_four_terminal import bounds
from terminal_commutators import parity, sign


@lru_cache(None)
def path_count(weight):
    if min(weight)<0 or sum(weight)%2 or max(weight[:2])>4:
        return 0
    paths = sorted(v for v in product(*(range(min(n,1)+1) if i<2 else range(n+1)
                                       for i,n in enumerate(weight))) if sum(v)%2)
    @lru_cache(None)
    def count(rem,first,slots):
        if slots==1:
            return int(sum(rem)%2 and max(rem[:2])<=1 and first<len(paths) and rem>=paths[first])
        total=0
        for i in range(first,len(paths)):
            v=paths[i]
            after=tuple(x-y for x,y in zip(rem,v))
            if min(after)>=0 and sum(after)>=slots-1 and max(after[:2])<=slots-1:
                total+=count(after,i+parity(v),slots-1)
        return total
    return count(weight,0,4)


@lru_cache(None)
def weyl_spectra():
    # det(1-z g) is the product of 1-s z^k over the signed permutation cycles.
    spectra = Counter()
    for perm in permutations(range(4)):
        for signs in product((-1,1),repeat=4):
            unused=set(range(4)); polynomial=[1]
            while unused:
                i=min(unused); length=0; orientation=1
                while i in unused:
                    unused.remove(i)
                    length+=1; orientation*=signs[i]; i=perm[i]
                poly=polynomial+[0]*length
                for j,n in enumerate(polynomial):
                    poly[j+length]-=orientation*n
                polynomial=poly
            spectra[(tuple(polynomial),prod(signs))]+=1
    assert sum(spectra.values())==384
    return spectra


@lru_cache(None)
def molien(weight,reflection):
    if min(weight)<0 or max(weight[:2])>4:
        return 0
    value=0
    for (poly,flip),size in weyl_spectra().items():
        symmetric=[1]
        for k in range(1,max(weight[2:])+1):
            symmetric.append(-sum(poly[j]*symmetric[k-j] for j in range(1,min(4,k)+1)))
        term=prod((-1)**n*poly[n] for n in weight[:2])*prod(symmetric[n] for n in weight[2:])
        value+=size*flip**reflection*term
    assert value%384==0 and value>=0
    return value//384


def mixed_pfaffian(matrices):
    value=0
    for mask in range(1,16):
        subset=[a for j,a in enumerate(matrices) if mask&(1<<j)]
        value+=(-1)**(4-len(subset))*pfaffian(sum(subset))
    return value


def absorption():
    rng=np.random.default_rng(230117)
    results=[]
    for _ in range(4):
        raw=rng.integers(-2,3,size=(5,8,8),dtype=np.int64)
        paths=[a-a.T for a in raw[:4]]
        B=raw[-1]
        lhs=int(np.trace(B))*mixed_pfaffian(paths)
        terms=[mixed_pfaffian(paths[:j]+[B@a+a@B.T]+paths[j+1:]) for j,a in enumerate(paths)]
        assert sum(terms)==lhs and all(terms)
        results.append(dict(left=lhs,pfaffian_terms=terms))
    return results


def run(eight,nine,fifth):
    mixed,scalar,pure=map(read,(eight,nine,fifth))
    assert mixed["classical_dimension"]==mixed["cartan_rank"]==3
    assert scalar["upper"]==scalar["lower"]==6 and int(scalar["minor"]["determinant"])!=0
    assert pure["boundary_and_augmented_ranks"]==[19,20]
    results=[]
    for t in range(6):
        w=(t,t)+(5-t,)*3
        even=odd=0
        shifts=[]
        for spin,ss in ((0,1),(1,-1)):
            for perm in permutations(range(3)):
                mu=(t+spin,t-spin)+tuple(5-t+perm[i]-i for i in range(3))
                coeff=ss*sign(perm)
                paths=path_count(mu); c=molien(mu,1)
                assert paths==c
                even+=coeff*molien(mu,0);odd+=coeff*c
                shifts.append(dict(weight=mu,sign=coeff,odd_paths=paths))
        formal=bounds(t)
        assert formal["cartan"]==even and odd>=0
        if t in (0,2,4):
            assert formal["formal"]==0
        results.append(dict(t=t,weight=w,length=sum(w),even_cartan=even,
                            odd_cartan=odd,odd_classical=odd,even_classical=even+int(t==5),
                            character_terms=shifts))
    return dict(group="SO(8)",level=30,degree_order="two gauginos, three scalars",
                pfaffian_absorption=absorption(),signed_permutations=384,sectors=results)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    base=Path(__file__).with_name("rank_four")
    p.add_argument("--eight",type=Path,default=base/"so8_mixed"/"so8_mixed_pairing.json")
    p.add_argument("--nine",type=Path,default=base/"scalar_terminal.json")
    p.add_argument("--fifth",type=Path,default=base/"so8_fifth.json")
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    assert __debug__
    cert=run(args.eight,args.nine,args.fifth)
    m.write(args.out,cert)
    print([(r["even_classical"],r["odd_classical"]) for r in cert["sectors"]])
