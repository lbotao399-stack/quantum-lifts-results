"""Integral spin and flavour singlet dimensions for the F4 fifth-power channels."""
import argparse
import json
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import factorial
from pathlib import Path

CARTAN = {
    "F4": ((2,-1,0,0),(-1,2,-2,0),(0,-1,2,-1),(0,0,-1,2)),
    "B4": ((2,-1,0,0),(-1,2,-1,0),(0,-1,2,-2),(0,0,-1,2)),
}
ZERO = (0,0,0,0)


def reflect(weight, root):
    i, row = root
    return tuple(x-weight[i]*y for x,y in zip(weight,row))


@lru_cache(None)
def roots(group):
    simple = tuple(enumerate(CARTAN[group]))
    known = set(CARTAN[group])
    pending = list(known)
    while pending:
        w = pending.pop()
        for s in simple:
            v = reflect(w,s)
            if v not in known:
                known.add(v)
                pending.append(v)
    assert len(known) == {"F4":48,"B4":32}[group]
    return tuple(sorted(known))


@lru_cache(None)
def chamber(group, weight):
    sign = 1
    while 0 not in weight:
        if min(weight) > 0:
            return sign, tuple(x-1 for x in weight)
        i = next(j for j,x in enumerate(weight) if x<0)
        weight = reflect(weight,(i,CARTAN[group][i]))
        sign = -sign
    return 0, ZERO


@lru_cache(None)
def multiply(group, highest, k):
    # Multiply a Weyl numerator by the k-th Adams operation on the adjoint.
    result = Counter({highest:4})
    for root in roots(group):
        sign, target = chamber(group, tuple(x+1+k*y for x,y in zip(highest,root)))
        if sign:
            result[target] += sign
    return tuple((h,c) for h,c in result.items() if c)


@lru_cache(None)
def decomposition(group, cycles):
    if not cycles:
        return {ZERO:1}
    out = Counter()
    for highest,c in decomposition(group,cycles[1:]).items():
        for target, n in multiply(group,highest,cycles[0]):
            out[target] += c*n
    return {h:c for h,c in out.items() if c}


@lru_cache(None)
def multiplicity(group, bosons, fermions):
    from charge_maps import content_patterns  # Loads the native character runtime on this path.
    from exact_dims import cycle_data
    den = 1
    cycles = {():1}
    for count, odd in [(n,False) for n in bosons]+[(n,True) for n in fermions]:
        den *= factorial(count)
        new = Counter()
        for a,c in cycles.items():
            for b,n in cycle_data(count,odd):
                new[tuple(sorted(a+b,reverse=True))] += c*n
        cycles = {u:n for u,n in new.items() if n}
    total = sum(c*decomposition(group,parts).get(ZERO,0) for parts,c in cycles.items())
    assert total >= 0 and total % den == 0
    return total//den


def singlets(group, mu, length):
    from charge_maps import content_patterns
    patterns = Counter()
    for ds, spin_sign in ((0,1),(1,-1)):
        for perm in permutations(range(3)):
            sign = spin_sign * (-1)**sum(i>j for n,i in enumerate(perm) for j in perm[n+1:])
            weight = (mu[0]+ds,mu[1]-ds)+tuple(mu[i+2]+perm[i]-i for i in range(3))
            for b,f,n in content_patterns(weight,length):
                patterns[b,f] += sign*n
    counts = Counter()
    for (b,f),n in patterns.items():
        if n:
            counts[sum(f)] += n*multiplicity(group,b,f)
    assert all(n>=0 for n in counts.values())
    return {d:n for d,n in sorted(counts.items()) if n}


def invariant_controls(group):
    degrees = (2,6,8,12) if group=="F4" else (2,4,6,8)
    max_degree = 12
    sym = [1]+[0]*max_degree
    ext = sym.copy()
    for degree in degrees:
        for k in range(degree,max_degree+1):
            sym[k] += sym[k-degree]
        for k in range(max_degree,2*degree-2,-1):
            ext[k] += ext[k-(2*degree-1)]
    for k in range(1,max_degree+1):
        assert multiplicity(group,(k,),()) == sym[k]
        assert multiplicity(group,(),(k,)) == ext[k]
    return dict(symmetric=sym, exterior=ext)


def run():
    from rank_four_s5 import dimensions
    comparisons = []
    for t in (5,4,3):
        mu = (t,t)+(5-t,)*3
        for length in range(13-t,16-t):
            actual = singlets("B4",mu,length)
            assert actual == dimensions("SO_odd",mu,length)
            comparisons.append(dict(weight=mu, length=length, dimensions=actual))
    controls = {g:invariant_controls(g) for g in CARTAN}
    rows = []
    for t in (5,4,3,2,1,0):
        mu = (t,t)+(5-t,)*3
        terminal = 15-t
        dimensions_by_length = {length:singlets("F4",mu,length) for length in range(2,terminal+1)}
        euler = sum((-1)**n*sum(v.values()) for n,v in dimensions_by_length.items())
        rows.append(dict(t=t,weight=mu,terminal_length=terminal,by_length=dimensions_by_length,
                         alternating_sum=euler))
        print(t,[sum(dimensions_by_length[n].values()) for n in (terminal-2,terminal-1,terminal)],
              euler,flush=True)
    return dict(group="F4", controls=controls, orthogonal_character_checks=comparisons, rows=rows)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert __debug__
    a.out.write_text(json.dumps(run(),indent=2)+"\n")
