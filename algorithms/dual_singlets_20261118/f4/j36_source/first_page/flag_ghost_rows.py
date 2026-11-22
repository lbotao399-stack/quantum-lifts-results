"""Reduce the 433 independent boundary directions modulo the zero-score ghost splittings."""
from functools import cache
from itertools import combinations, product
import json
from random import Random

import numpy as np
from sympy import Matrix

from flag_boundaries import HERE, polynomial, recipe, trace_words
from ghost_trace_quotient import rotate, sign


@cache
def necklaces(word):
    f, b = rotate(word), rotate(word[::-1])
    phase = (-1)**len(word)
    if f == b and phase < 0:
        return None,0
    return min(f,b), (1 if f <= b else phase)


def moments(count):
    words = sorted({necklaces(w)[0] for n in range(7) for w in product(range(3),repeat=n)
                    if necklaces(w)[1]})
    rng = Random(644331)
    samples = {w:[rng.randrange(-2,3) for _ in range(count)] for w in words}
    # The final column commutes the scalars as well. Every marked commutator
    # must then vanish, with the same reversal convention as in the ring.
    for w in words:
        val = samples[tuple(sorted(w))][0] if len(w)%2 == 0 else 0
        samples[w] = np.array([*samples[w],val],dtype=np.int64)
    samples[()] = np.ones(count+1,dtype=np.int64)
    return samples


def catalogue_values(samples):
    _, reps, size = trace_words()
    ncols = len(samples[()])
    result = np.zeros((size,4,4,ncols),dtype=np.int16)
    offset = 1
    for length, codes in reps.items():
        for ix, code in enumerate(codes):
            word = tuple(int(code >> (2*i)) & 3 for i in reversed(range(length)))
            for a,b in combinations(range(length),2):
                rest = tuple(x for p,x in enumerate(word) if p not in (a,b))
                if 3 in rest:
                    continue
                key,z = necklaces(rest)
                if z:
                    i,j = word[a],word[b]
                    result[offset+ix,i,j] += z*samples[key]
                    result[offset+ix,j,i] += z*samples[key]
        offset += len(codes)
    return result


PAIRS = tuple(combinations(range(4),2))


def paired(values,left,right):
    # Antisymmetrize the two row slots and the two column slots separately.
    # Each half is shared by many monomials with two short traces.
    out = {}
    for i,j in PAIRS:
        for a,b in PAIRS:
            out[(i,j),(a,b)] = (values[left,i,a].astype(np.int32)*values[right,j,b]
                                 -values[left,i,b].astype(np.int32)*values[right,j,a]
                                 -values[left,j,a].astype(np.int32)*values[right,i,b]
                                 +values[left,j,b].astype(np.int32)*values[right,i,a])
    return out


def row(seed,values):
    coefficients,ids = polynomial(seed)
    assert ids.shape[0] == 4
    left, li = np.unique(ids[:2].T,axis=0,return_inverse=True)
    right, ri = np.unique(ids[2:].T,axis=0,return_inverse=True)
    first, second = paired(values,*left.T), paired(values,*right.T)
    out = np.zeros(values.shape[-1],dtype=np.int64)
    for top in PAIRS:
        bottom = tuple(i for i in range(4) if i not in top)
        for a in PAIRS:
            b = tuple(i for i in range(4) if i not in a)
            phase = sign(top+bottom)*sign(a+b)
            scores = first[top,a][li].astype(np.int64)*second[bottom,b][ri]
            out += phase*(coefficients @ scores)
    return out.tolist()


def calculate():
    saved = json.loads((HERE/"flag_433.json").read_bytes())
    rows = saved["rows"]
    active = [(i,s) for i,(family,s) in enumerate(rows) if family=="short" and len(recipe(s)[0])==4]
    values = catalogue_values(moments(8))
    matrix = [row(seed,values) for _,seed in active]
    assert all(r[-1] == 0 for r in matrix)
    data = Matrix([r[:-1] for r in matrix])
    chosen = list(data.T.rref()[1])
    points = list(data[chosen,:].rref()[1])
    square = data.extract(chosen,points)
    assert len(chosen) == len(points) == 4 and square.det()
    coords = data.extract(range(len(active)),points)*square.inv()
    assert coords*data[chosen,:] == data
    return dict(number_of_boundary_rows=len(rows),four_trace_rows=len(active),
                moment_basis_rows=[active[i][0] for i in chosen],moment_basis_points=points,
                determinant=str(square.det()),minor=[[int(x) for x in r] for r in square.tolist()],
                four_trace_recipes=[dict(index=i,recipe=s,coordinates=list(map(str,coords.row(j))))
                                    for j,(i,s) in enumerate(active)],
                dimension_with_explicit_ghost_preimages=len(rows)-len(chosen),
                upper_bound_on_restricted_dual_unknowns=704-len(rows)+len(chosen))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
