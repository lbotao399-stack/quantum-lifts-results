"""Triality on the lifted fermion modes of a charge-one SO(8) instanton."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
import json


def subset_weights(weights, degree):
    return Counter(tuple(map(sum, zip(*s))) if s else (0,)*len(weights[0])
                   for s in combinations(weights, degree))


def highest(character, rank):
    out = {}
    for weight in character:
        if min(weight) < 0:
            continue
        coefficient = sum((-1)**sum(s)*character.get(tuple(w+2*d for w,d in zip(weight,s)),0)
                          for s in product(range(2), repeat=rank))
        if coefficient:
            out[weight] = coefficient
    return out


@lru_cache(None)
def exterior(degree, rotated):
    if not rotated:
        weights = list(product((-1,1), repeat=3))
        ans = highest(subset_weights(weights, degree),3)
        assert all(n>0 for n in ans.values())
        assert sum(n*prod(i+1 for i in w) for w,n in ans.items()) == comb(8,degree)
        return ans
    # The cyclic permutation has two fixed weights and two three-cycles.
    factors = ((1,1), (1,-1), (3,1), (3,-1))
    ch = Counter()
    for k in range(5):
        for subset in combinations(factors,k):
            if sum(d for d,w in subset)==degree:
                ch[(sum(w for d,w in subset),)] += 1
    # Weyl differences can be negative for a twined character.
    return {(i,):n for i in range(5)
            if (n:=ch[(i,)]-ch[(i+2,)])}


@lru_cache(None)
def tensor_degrees(first, second, rotated):
    out = Counter()
    for a,v in exterior(first,rotated).items():
        for b,w in exterior(second,rotated).items():
            factors = [range(abs(i-j),i+j+1,2) for i,j in zip(a,b)]
            for irrep in product(*factors):
                out[irrep] += v*w
    return out


@lru_cache(None)
def invariant_trace(degrees, rotated):
    a,b,c,d = sorted(degrees)
    left = tensor_degrees(a,b,rotated)
    right = tensor_degrees(c,d,rotated)
    return sum(n*right.get(w,0) for w,n in left.items())


def standard(degrees):
    identity = invariant_trace(tuple(degrees), False)
    cycle = invariant_trace(tuple(degrees), True)
    multiplicity, remainder = divmod(identity-cycle, 3)
    assert not remainder and multiplicity>=0
    return multiplicity


def run():
    # Each adjoint Weyl fermion has four exact and eight lifted Grassmann variables.
    # The two insertions have eight scalars and two matter fermions of opposite chirality.
    rows=[]
    for exact in product(range(5), repeat=3):
        left=(6,)+tuple(4-n for n in exact)
        right=(0,)+tuple(2+n for n in exact)
        assert all(a+b==6 for a,b in zip(left,right))
        l,r=standard(left),standard(right)
        if l and r:
            rows.append(dict(left_exact=[4,*exact],left_lifted=left,right_lifted=right,
                             left_standard=l,right_standard=r))
    first={d:sum(standard(v) for v in product(range(d+1),repeat=4) if sum(v)==d)
           for d in range(5)}
    assert first=={0:0,1:0,2:0,3:0,4:20}
    counts=Counter((sum(x["left_lifted"]),sum(x["right_lifted"])) for x in rows)
    return dict(exact_modes=16,lifted_modes=32,leading_insertions=[[10,4,4,4],[0,6,6,6]],
                quartic_insertions=2,first_standard_multiplicities=first,
                allowed_assignments=len(rows),by_lifted_degree=[[*k,n] for k,n in sorted(counts.items())],
                assignments=rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args=parser.parse_args()
    result=run()
    args.out.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(result["allowed_assignments"],result["by_lifted_degree"])
