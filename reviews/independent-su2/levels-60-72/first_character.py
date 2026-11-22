"""Exact SU(2) colour and SU(3) flavour singlets in the spin-zero first target."""
from argparse import ArgumentParser
from collections import Counter
from itertools import combinations_with_replacement, permutations
import json
from math import comb

import numpy as np


def calculate(n):
    m=3+2*n
    limit=7+2*n+1
    counts=np.zeros((m+3,m+3,m+3,limit+1,2*limit+1),dtype=np.int64)
    counts[0,0,0,0,limit]=1
    assert 2**9*comb(limit+12,12)<2**62
    letters=[(tuple((s>>i)&1 for i in range(3)),z) for s in range(8) for z in (-1,0,1)]
    for w,z in letters:
        if not any(w):
            continue
        power=1
        while power<=limit and all(power*i<=m+2 for i in w):
            shift=tuple(power*i for i in (*w,1,z))
            dst=tuple(slice(max(i,0),size+min(i,0)) for size,i in zip(counts.shape,shift))
            src=tuple(slice(max(-i,0),size-max(i,0)) for size,i in zip(counts.shape,shift))
            counts[dst]+=counts[src].copy()
            if sum(w)==2:
                break
            power*=2
    pairs=Counter()
    for (w,z),(v,u) in combinations_with_replacement(letters,2):
        if (w,z)==(v,u) and sum(w)%2:
            continue
        pairs[tuple(i+j for i,j in zip(w,v)),z+u]+=1
    result={}
    for y in range(7+2*n-1,limit+1):
        complete=[0,0]
        for perm in permutations(range(3)):
            sign=(-1)**sum(i>j for k,i in enumerate(perm) for j in perm[k+1:])
            for (w,z),multiplicity in pairs.items():
                powers=tuple(m+p-i-w[i] for i,p in enumerate(perm))
                if min(powers)<0:
                    continue
                for colour in (0,1):
                    complete[colour]+=sign*multiplicity*int(counts[(*powers,y-2,limit+colour-z)])
        result[y]=dict(flavour_singlet_at_colour_weights=complete,dimension=complete[0]-complete[1])
    return result


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('n',nargs='+',type=int)
    args=parser.parse_args()
    if not __debug__ or min(args.n)<1:
        parser.error('enable assertions and use positive n')
    print(json.dumps({n:calculate(n) for n in args.n},indent=2))
