"""Exact BMN cochain dimensions and Euler indices after removing the three raising supercharges."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import permutations, product
import json
from math import factorial
from pathlib import Path

from b2_schur_complex import centralizer, colour_trace, partitions


@lru_cache(None)
def g2_dominant(weight):
    first, second = weight
    while min(first, second) < 0:
        if first < 0:
            first, second = -first, second + first
        else:
            first, second = first + 3*second, -second
    return first, second


def g2_reflections(w):
    a, b = w
    return (-a, b+a), (a+3*b, -b)


def g2_orbit(start):
    signs = {start: 1}
    pending = [start]
    while pending:
        w = pending.pop()
        for v in g2_reflections(w):
            if v not in signs:
                signs[v] = -signs[w]
                pending.append(v)
    return signs


G2_ROOTS = tuple(sorted(set(g2_orbit((2,-1))) | set(g2_orbit((-3,2)))))
G2_WEYL = g2_orbit((1,1))


@lru_cache(None)
def g2_weight(cycles, weight):
    if not cycles:
        return int(weight == (0,0))
    k = cycles[-1]
    a, b = weight
    if a*2+b*3 > 3*sum(cycles) or a+2*b > 2*sum(cycles):
        return 0
    tail = cycles[:-1]
    return 2*g2_weight(tail,weight) + sum(
        g2_weight(tail,g2_dominant((a-k*i,b-k*j))) for i,j in G2_ROOTS)


def g2_trace(cycles):
    return sum(sign*g2_weight(cycles,g2_dominant((1-a,1-b)))
               for (a,b),sign in G2_WEYL.items())


@lru_cache(None)
def su2_colour_weight(cycles, weight):
    """Adams character of the adjoint, in units of the SU(2) root."""
    if not cycles:
        return int(weight == 0)
    if weight > sum(cycles):
        return 0
    k = cycles[-1]
    return sum(su2_colour_weight(cycles[:-1], abs(weight-shift)) for shift in (-k, 0, k))


def colour(cycles, group):
    if group == "SU(2)":
        return (su2_colour_weight(cycles, 0) - su2_colour_weight(cycles, 1),)
    if group == "USp(4)":
        return (colour_trace(cycles),)
    if group == "G2":
        return (g2_trace(cycles),)
    if group in ("F4", "B4"):
        from f4_singlet_characters import decomposition, ZERO
        return (decomposition(group,cycles).get(ZERO,0),)
    if group == "E6":
        from bmn_primary.exceptional_trace import e6_traces
        return e6_traces(cycles)
    if group in ("E7","E8"):
        from bmn_primary.exceptional_trace import invariant_trace
        return (invariant_trace(group,cycles),)
    if group != "SU(3)":
        raise ValueError(f"Unsupported gauge group {group}")
    from unitary_characters import adjoint_traces
    trace, twisted = adjoint_traces(tuple(reversed(cycles)), 3)
    if (trace-twisted) % 2:
        raise ArithmeticError("Nonintegral traces on the two colour eigenspaces")
    return ((trace+twisted)//2, (trace-twisted)//2)


def flavour_characters(maximum):
    @lru_cache(None)
    def polynomial(cycles):
        if not cycles:
            return {(0,0,0): 1}
        k, tail = cycles[-1], cycles[:-1]
        odd = (-1)**(k-1)
        letters = ((k,0,0,1),(0,k,0,1),(0,0,k,1),
                   (k,k,0,odd),(0,k,k,odd),(k,0,k,odd),(k,k,k,1))
        terms = defaultdict(int)
        for (a,b,c), coefficient in polynomial(tail).items():
            for i,j,h,sign in letters:
                if a+b+c+i+j+h <= maximum:
                    terms[a+i,b+j,c+h] += sign*coefficient
        return {w:v for w,v in terms.items() if v}
    return polynomial


def cochains(maximum, group):
    poly = flavour_characters(maximum)
    rho = (2,1,0)
    alternant = [(tuple(a-b for a,b in zip(rho,w)),
                  (-1)**sum(a<b for i,a in enumerate(w) for b in w[i+1:])) for w in permutations(rho)]
    shapes = [(a,b,c) for a in range(maximum+1) for b in range(min(a,maximum-a)+1)
              for c in range(min(b,maximum-a-b)+1)]
    rows = {}
    for length in range(maximum+1):
        scale = factorial(length)
        accumulators = [defaultdict(int) for _ in range(2 if group in ("SU(3)","E6") else 1)]
        for cycles in partitions(length):
            coefficients = colour(cycles,group)
            for acc, multiplicity in zip(accumulators,coefficients):
                if not multiplicity:
                    continue
                factor = (scale//centralizer(cycles))*multiplicity
                for w,n in poly(cycles).items():
                    if w[0]>=w[1]>=w[2]:
                        acc[w] += factor*n
        result = {}
        for shape in shapes:
            if max(shape)>length or sum(shape)<length or sum(shape)>3*length:
                continue
            values = []
            for acc in accumulators:
                total = 0
                for shift, sign in alternant:
                    w = tuple(a+b for a,b in zip(shape,shift))
                    if min(w)>=0:
                        total += sign*acc.get(tuple(sorted(w,reverse=True)),0)
                value,remainder = divmod(total,scale)
                if value<0 or remainder:
                    raise ArithmeticError((length,shape,total,scale))
                values.append(value)
            if any(values):
                result[shape] = values
        rows[length] = result
    return shapes, rows


def primary_indices(maximum, group):
    shapes, levels = cochains(maximum,group)
    parities = 2 if group in ("SU(3)","E6") else 1
    rows = []
    for shape in shapes:
        a,b,c = shape
        if c<2:
            continue
        strips = [(w,(-1)**(a+b+c-sum(w))) for w in product(range(b,a+1),range(c,b+1),range(c+1))]
        dimensions = {}
        ordinary = [0]*parities
        index = [0]*parities
        for length, chars in levels.items():
            for k,v in enumerate(chars.get(shape,())):
                ordinary[k] += (-1)**length*v
            count = [0]*parities
            for w,sign in strips:
                for k,v in enumerate(chars.get(w,())):
                    count[k] += sign*v
            if min(count)<0:
                raise ArithmeticError(("negative lowest-weight multiplicity",length,shape,count))
            if any(count):
                dimensions[length] = count
                for k,v in enumerate(count):
                    index[k] += (-1)**length*v
        rows.append(dict(weight=shape,level=2*sum(shape),cochains_by_length=dimensions,
                         euler=ordinary,primary_euler=index))
    return dict(group=group,maximum_level=2*maximum,colour_components=(["even","odd"] if parities==2 else ["all"]),
                rows=sorted(rows,key=lambda r:(r["level"],r["weight"])))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--group",choices=("SU(2)","USp(4)","G2","SU(3)","F4","B4","E6","E7","E8"),required=True)
    p.add_argument("--maximum-level",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    args = p.parse_args()
    if args.maximum_level<0 or args.maximum_level%2:
        p.error("Choose a nonnegative even BMN level")
    if len(G2_ROOTS)!=12 or len(G2_WEYL)!=12 or any(
            g2_trace(cycles)!=expected for cycles,expected in [((),1),((1,),0),((1,1),1),((2,),1)]):
        raise ArithmeticError("Inconsistent G2 root system")
    data = primary_indices(args.maximum_level//2,args.group)
    args.out.write_text(json.dumps(data,separators=(",",":"))+"\n")
