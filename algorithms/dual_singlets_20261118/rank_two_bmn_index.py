"""Find odd-length BMN classes using exact colour characters and protected highest weights."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import random

from flint import nmod_mat
import mixed_incoming as m
from sp3_omega_products import all_products, canonical, cartan, letter_odd, slot_odd
from trace_descent import dimensions


def colour_profiles(weight, length):
    degree = sum(weight)-length
    profiles = set()
    for strengths in range(degree//2+1):
        fermions = degree-2*strengths
        scalars = length-fermions-strengths
        for n in range(scalars+1):
            for k in range(max(0,fermions-3),min(3,fermions)+1):
                for f in range(strengths+1):
                    profile = (n,k,f), (scalars-n,fermions-k,strengths-f)
                    if all(sum(p)%2 == 0 for p in profile):
                        profiles.add(tuple(sorted(profile)))
    return sorted(profiles)


@lru_cache(None)
def slot_words(profile):
    n, k, f = profile
    result = defaultdict(list)
    for a in range(n+1):
        for b in range(n-a+1):
            for subset in combinations((3,5,6), k):
                masks = (1,)*a+(2,)*b+(4,)*(n-a-b)+subset+(7,)*f
                charge = tuple(sum((mask>>i)&1 for mask in masks) for i in range(3))
                result[charge].append(tuple((0,0,v) for v in sorted(masks)))
    return result


@lru_cache(maxsize=350)
def monomials(weight, profiles):
    left, right = map(slot_words, profiles)
    keys = set()
    for w, words in left.items():
        complement = tuple(a-b for a,b in zip(weight,w))
        for x in words:
            for y in right.get(complement, ()):
                s, v = canonical((x,y), slot_odd)
                if s:
                    keys.add(v)
    return sorted(keys)


def lower(monomial, axis):
    result = Counter()
    for i, slot in enumerate(monomial):
        for j, (_, _, mask) in enumerate(slot):
            if not (mask & (1<<axis)) or mask & (2<<axis):
                continue
            jet = (0,0,mask^(3<<axis))
            s, v = canonical(slot[:j]+(jet,)+slot[j+1:], letter_odd)
            if s:
                t, x = canonical(monomial[:i]+(v,)+monomial[i+1:], slot_odd)
                if t:
                    result[x] += s*t
    return {k:v for k,v in result.items() if v}


def multiplicity(weight, profiles):
    rho = (2,1,0)
    value = 0
    for p in permutations(rho):
        sign = (-1)**sum(x<y for i,x in enumerate(p) for y in p[i+1:])
        shifted = tuple(n+a-b for n,a,b in zip(weight,rho,p))
        value += sign*len(monomials(shifted,profiles))
    return value


def quotient(weight,length):
    # V_weight / (F1 V_(weight+alpha1) + F2 V_(weight+alpha2)).
    # A character multiplicity fixes the rational denominator rank in every block.
    maps, seen, quotdim = {}, 0, 0
    for profiles in colour_profiles(weight,length):
        words = monomials(weight,profiles)
        if not words:
            continue
        size = multiplicity(weight,profiles)
        assert 0 <= size <= len(words)
        seen += len(words)
        if not size:
            continue
        index = {w:j for j,w in enumerate(words)}
        rows = []
        for axis in range(2):
            above = list(weight)
            above[axis] += 1
            above[axis+1] -= 1
            for w in monomials(tuple(above),profiles):
                row = [0]*len(words)
                for k,n in lower(w,axis).items():
                    row[index[k]] = n
                rows.append(row)
        field = nmod_mat(rows,m.P) if rows else nmod_mat(0,len(words),m.P)
        reduced, rank = field.rref()
        assert rank == len(words)-size
        piv, j = [], 0
        for i in range(rank):
            while not reduced[i,j]:
                j += 1
            piv.append(j)
            j += 1
        free = sorted(set(range(len(words)))-set(piv))
        for i,q in enumerate(free):
            maps[words[q]] = [(quotdim+i,1)]
        for i,j in enumerate(piv):
            row = [(quotdim+k, -int(reduced[i,q]) % m.P) for k,q in enumerate(free) if reduced[i,q]]
            if row:
                maps[words[j]] = row
        quotdim += size
    return maps,seen,quotdim


def highest_products(weight,length,sample_count,seed):
    maps, weight_dim, highest = quotient(weight,length)
    piv, witnesses = {}, []
    if highest:
        families = all_products(4,(0,0,*weight),length)
        rows = [w for group in families.values() for w in group]
        random.Random(seed).shuffle(rows)
        for w in rows[:sample_count]:
            image = Counter()
            for monomial,n in cartan(w,2).items():
                for j,c in maps.get(monomial,()):
                    image[j] += c*n
            reduced = {j:c%m.P for j,c in image.items() if c%m.P}
            while reduced:
                lead = min(reduced)
                if lead not in piv:
                    inv = pow(reduced[lead],-1,m.P)
                    piv[lead] = {j:c*inv%m.P for j,c in reduced.items()}
                    witnesses.append(w)
                    break
                multiple = reduced[lead]
                for j,c in piv[lead].items():
                    value = (reduced.get(j,0)-multiple*c)%m.P
                    if value:
                        reduced[j] = value
                    else:
                        reduced.pop(j,None)
            if len(piv) == highest:
                break
    return dict(length=length,cartan_weight_dimension=weight_dim,cartan_highest_dimension=highest,
                product_lower_bound=len(piv),products=witnesses,pivots=sorted(piv))


def run(weight,sample_count,output):
    weight = tuple(weight)
    assert min(weight)>=0 and tuple(sorted(weight,reverse=True))==weight
    centre = (0,0,*weight)
    dims = {n:sum(dimensions(centre,n,4,"Sp").values()) for n in range(max(weight),sum(weight)+1)}
    euler = sum((-1)**n*d for n,d in dims.items())
    records=[]
    for n,dim in dims.items():
        if n%2 or not dim:
            continue
        records.append(highest_products(weight,n,sample_count,924917+127*n))
        print(weight,n,records[-1]["product_lower_bound"],records[-1]["cartan_highest_dimension"],flush=True)
    count = sum(row["product_lower_bound"] for row in records)
    result = dict(group="USp(4)",weight=weight,level=2*sum(weight),prime=m.P,
                  dimensions=dims,euler_with_length_sign=euler,even_product_lower_bound=count,
                  odd_cohomology_lower_bound=max(0,count-euler),sample_per_length=sample_count,
                  even_length_records=records)
    m.write(output,result)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight",nargs=3,type=int,required=True)
    parser.add_argument("--samples",type=int,default=256)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    record=run(args.weight,args.samples,args.out)
    print("Odd cohomology lower bound",record["odd_cohomology_lower_bound"],flush=True)
