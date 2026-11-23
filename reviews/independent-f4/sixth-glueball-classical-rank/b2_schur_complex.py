"""Exact USp(4) cochain multiplicities of polynomial Taylor representations."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
import json
from math import factorial, prod
from pathlib import Path


@lru_cache(None)
def partitions(n, minimum=1):
    if not n:
        return ((),)
    return tuple((k,)+tail for k in range(minimum,n+1)
                 for tail in partitions(n-k,k))


@lru_cache(None)
def centralizer(cycles):
    return prod(k**m * factorial(m) for k,m in Counter(cycles).items())


def transpose(shape):
    return tuple(sum(n >= j for n in shape) for j in range(1,max(shape,default=0)+1))


@lru_cache(None)
def dimension(shape):
    dual = transpose(shape)
    hooks = prod(n-j+dual[j]-i-1 for i,n in enumerate(shape) for j in range(n))
    return factorial(sum(shape))//hooks


@lru_cache(None)
def character(shape, cycles):
    """Murnaghan–Nakayama, moving one bead down by the longest cycle."""
    if not cycles or cycles[-1] == 1:
        return dimension(shape)
    length, step = len(shape), cycles[-1]
    beads = tuple(n+length-i-1 for i,n in enumerate(shape))
    result = 0
    for i, bead in enumerate(beads):
        destination = bead-step
        if destination < 0 or destination in beads:
            continue
        height = sum(destination < other < bead for other in beads)
        moved = sorted(beads[:i]+beads[i+1:]+(destination,), reverse=True)
        smaller = tuple(b-length+j+1 for j,b in enumerate(moved) if b-length+j+1)
        result += (-1)**height * character(smaller,cycles[:-1])
    return result


C2_ROOTS = ((2,0),(-2,0),(0,2),(0,-2),(1,1),(1,-1),(-1,1),(-1,-1))


@lru_cache(None)
def colour_weight(cycles, u, v):
    if u+v > 2*sum(cycles) or (u+v)%2:
        return 0
    if not cycles:
        return int(u==v==0)
    k, tail = cycles[-1], cycles[:-1]
    count = 2*colour_weight(tail,u,v)
    for a,b in C2_ROOTS:
        weight = sorted((abs(u-k*a), abs(v-k*b)), reverse=True)
        count += colour_weight(tail,*weight)
    return count


@lru_cache(None)
def colour_trace(cycles):
    """Trace of a permutation on the invariant adjoint tensor power."""
    total = 0
    for perm in permutations((2,1)):
        for u,v in product((-1,1), repeat=2):
            weight = sorted((abs(2-u*perm[0]),abs(1-v*perm[1])),reverse=True)
            total += (1 if perm==(2,1) else -1)*u*v*colour_weight(cycles,*weight)
    return total


@lru_cache(maxsize=12000)
def jet_character(cycles, excess):
    """Power sums for positive symmetric jets along each cycle of the letters."""
    if not excess:
        return ((cycles, Fraction(1)),)
    if not cycles:
        return ()
    k, tail = cycles[-1], cycles[:-1]
    values = Counter()
    for extra in range(excess//k+1):
        for small in partitions(1+extra):
            lifted = tuple(k*i for i in small)
            reciprocal = Fraction(1,centralizer(small))
            for rest, c in jet_character(tail,excess-k*extra):
                values[tuple(sorted(lifted+rest))] += c*reciprocal
    return tuple(values.items())


@lru_cache(None)
def multiplicity(shape, length, trace=colour_trace):
    excess = sum(shape)-length
    if excess < 0:
        return 0
    ans = Fraction(0)
    for cycles in partitions(length):
        colour = trace(cycles)
        if not colour:
            continue
        pairing = sum(c*character(shape,jet_cycles)
                      for jet_cycles,c in jet_character(cycles,excess))
        ans += Fraction((-1)**(length-len(cycles))*colour,centralizer(cycles))*pairing
    assert ans.denominator==1 and ans>=0
    return int(ans)


def partner(weight, order):
    # All three odd-coordinate degrees in this calculation are at least two.
    # The first two rows of the source diagram both have three boxes.
    assert len(weight)==3 and tuple(sorted(weight,reverse=True))==weight and min(weight)>=2
    tail = tuple(k-2-order for k in weight)
    if min(tail)<0:
        return None
    return (3+order,3+order)+transpose(tail)


def calculate():
    here = Path(__file__).resolve().parent/"b2_level36"
    controls = []
    records = []
    for filename in ("rational_products.json", "products_10_5_5.json", "products_11_5_5.json"):
        data = json.loads((here/filename).read_bytes())
        weight = tuple(data["weight"])
        dims = {int(k):sum(r.values()) for k,r in data["cochains_by_length_and_fermion_degree"].items()}
        shape = partner(weight,0)
        assert shape==transpose(weight)
        computed = {y:multiplicity(shape,y) for y in dims}
        assert computed==dims
        controls.append(dict(weight=weight,shape=shape,cochains=computed))
        for order in (1,2,3,4):
            shape = partner(weight,order)
            if shape is None:
                records.append(dict(weight=weight,order=order,shape=None,cochains={}))
                continue
            sizes = {y:multiplicity(shape,y) for y in range(max(2,min(dims)-order),max(dims)+2-order)}
            records.append(dict(weight=weight,order=order,shape=shape,cochains=sizes))
    return dict(group="USp(4)", colour_roots=C2_ROOTS, bmn_character_controls=controls,
                order_selection="filtered superconformal equivariance",
                partner_multiplicities=records)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    args.out.write_text(json.dumps(calculate(),indent=2)+"\n")
