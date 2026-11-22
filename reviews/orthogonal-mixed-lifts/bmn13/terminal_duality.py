"""Cartan multiplicities at matched physical charges of SO(2n+1) and USp(2n)."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations,product
import json
from math import factorial,prod
from pathlib import Path


def partitions(n,low=1):
    if n==0:
        yield ()
    for j in range(low,n+1):
        for tail in partitions(n-j,j):
            yield (j,)+tail


@lru_cache(None)
def species(cycles,k,odd):
    coef=[1]+[0]*k
    for m,s in cycles:
        if odd:
            for j in range(k,m-1,-1):
                coef[j]-=(-1)**m*s*coef[j-m]
        else:
            for j in range(m,k+1):
                coef[j]+=s*coef[j-m]
    return coef[k]


@lru_cache(None)
def molien(n,weight):
    if min(weight)<0:
        return 0
    total=Fraction()
    for partition in partitions(n):
        z=prod(m**j*factorial(j) for m,j in Counter(partition).items())
        for signs in product((-1,1),repeat=len(partition)):
            cycles=tuple(sorted(zip(partition,signs)))
            total+=Fraction(prod(species(cycles,k,i<2) for i,k in enumerate(weight)),2**len(partition)*z)
    assert total.denominator==1 and total>=0
    return int(total)


def orbit_count(n,weight):
    import numpy as np
    if min(weight)<0:
        return 0
    shape=tuple(k+1 for k in (n,)+weight)
    counts=np.zeros(shape,dtype=np.int64)
    counts[(0,)*6]=1
    for g in product(*(range(min(k,1)+1) if i<2 else range(k+1) for i,k in enumerate(weight))):
        if not any(g) or sum(g)%2:
            continue
        lower=tuple(slice(0,k+1-v) for k,v in zip(weight,g))
        upper=tuple(slice(v,k+1) for k,v in zip(weight,g))
        # Identical odd site types have zero orbit sum. Even sites can repeat.
        for sites in (range(n,0,-1) if (g[0]+g[1])%2 else range(1,n+1)):
            counts[(sites,)+upper]+=counts[(sites-1,)+lower]
    return int(sum(counts[(k,)+weight] for k in range(n+1)))


def highest(n,w):
    total=0
    for spin in (0,1):
        for perm in permutations(range(3)):
            sign=(-1)**(spin+sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3)))
            shifted=(w[0]+spin,w[1]-spin)+tuple(w[2+i]+perm[i]-i for i in range(3))
            total+=sign*molien(n,shifted)
    return total


def run(n,extra):
    rows=[]
    for h in range(n+1):
        w=(h,h,n+extra-h,n-h,n-h)
        count=molien(n,w)
        assert count==orbit_count(n,w)
        rows.append(dict(degrees=w,length=sum(w),weight_multiplicity=count,highest_multiplicity=highest(n,w)))
    assert all(r['highest_multiplicity']>=0 for r in rows)
    return dict(rank=n,doubled_physical_charges=[0,0,2*n+2*extra,2*n,2*n],
                rows=rows,all_weights_total=sum(r['weight_multiplicity'] for r in rows),
                highest_total=sum(r['highest_multiplicity'] for r in rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rank',type=int,required=True)
    p.add_argument('--extra',type=int,default=2)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    assert __debug__ and 2<=a.rank<=6 and 0<=a.extra<=4
    result=run(a.rank,a.extra)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print([(r['weight_multiplicity'],r['highest_multiplicity']) for r in result['rows']],
          result['all_weights_total'],result['highest_total'])
