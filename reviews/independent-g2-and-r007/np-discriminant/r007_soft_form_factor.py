"""Four-scalar massless form factor of Dot's R-007 target at a regular SU(3) vacuum."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
import json
from math import factorial, prod
from pathlib import Path
from r007_free_radial import target


def four_scalar(words,vev,B,C):
    # Fourier factors i give i^4=1. All external lambda spinors have lambda_+=1.
    # The two B legs have only minus dotted momentum, the two C legs only plus.
    colours = (vev,B,C)
    polynomial = Counter()
    word_census = Counter()
    for w,c in words.items():
        legs = [l for tr in w for l in tr]
        if any(len(h)!=1 for _,_,h in legs):
            continue
        flavours = [h[0] for _,_,h in legs]
        if tuple(flavours.count(j) for j in range(3)) != (5,2,2):
            continue
        if any((fl==0 and (p or q)) or (fl==1 and p) or (fl==2 and q)
               for (p,q,h),fl in zip(legs,flavours)):
            continue
        colour=prod(sum(prod(colours[h[0]][k] for _,_,h in tr) for k in range(3)) for tr in w)
        derivatives=[[p+q for (p,q,h),fl in zip(legs,flavours) if fl==j] for j in (1,2)]
        factor=c*colour/prod(factorial(n) for slots in derivatives for n in slots)
        for b in permutations(derivatives[0]):
            for d in permutations(derivatives[1]):
                assert sum(b+d)==4
                polynomial[b+d]+=factor
        shape=tuple(tuple(sorted(x)) for x in derivatives)
        word_census[shape]+=1
    return {k:n for k,n in polynomial.items() if n},word_census


def trace_of(*matrices):
    return sum(prod(vals) for vals in zip(*matrices))


def gram_coordinate(vev,B,C):
    return (4*Fraction(56,15)*trace_of(vev,vev,vev)*trace_of(vev,vev)
            *(trace_of(B,B)*trace_of(C,C)-trace_of(B,C)**2))


def run():
    words=target()
    A,B,C=(1,2,-3),(1,0,-1),(0,1,-1)
    full,census=four_scalar(words,A,B,C)
    assert full[1,1,1,1]==gram_coordinate(A,B,C)==Fraction(-56448,5)
    assert len(set(A))==3 and sum(A)==sum(B)==sum(C)==0
    assert prod((A[i]-A[j])**2 for i in range(3) for j in range(i+1,3))==400
    cases=[]
    for a,b,c in [(A,B,C),((1,1,-2),B,C),(A,B,B),((-2,3,-1),B,C),
                   ((0,2,-2),B,C),(A,(1,-1,0),(-2,1,1))]:
        row,_=four_scalar(words,a,b,c)
        z=row.get((1,1,1,1),0)
        assert z==gram_coordinate(a,b,c)
        cases.append(dict(vev=a,first=b,second=c,coefficient=str(z)))
    # Each single-particle SU(4) weight is a scalar, a gaugino or a gauge field.
    scalar=[tuple(2*sign*int(i==j) for i in range(3)) for j in range(3) for sign in (1,-1)]
    gauginos=list(product((-1,1),repeat=3))
    weights=scalar+gauginos+[(0,0,0)]*2
    assert len(weights)==16 and max((w[1]+w[2]) for w in weights)==2
    assert not any(sum(w[1]+w[2] for w in state)==8
                   for n in range(4) for state in product(weights,repeat=n))
    return dict(terms=len(words),doubled_R=[10,4,4],doubled_left_spin=4,vev_discriminant=400,
                polynomial=[dict(degrees=k,coefficient=str(n)) for k,n in sorted(full.items())],
                scalar_trace_word_counts=[dict(derivatives=k,count=n) for k,n in sorted(census.items())],
                controls=cases,max_external_R2_plus_R3=1,minimum_target_external_particles=4)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        raise RuntimeError('enable assertions')
    r=run()
    a.out.write_text(json.dumps(r,indent=2)+'\n')
    print(r['polynomial'],flush=True)
