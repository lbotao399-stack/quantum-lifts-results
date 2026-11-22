"""Ghost placement on closed F4 colour graphs and a native six-vertex test."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
import gzip
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from math import gcd, lcm
from pathlib import Path
import numpy as np
from sympy import Matrix

from single_trace_flag import native26

HERE = Path(__file__).resolve().parent


def incidence(vertices, edges):
    matrix = np.zeros((vertices,len(edges)), dtype=int)
    for i,(u,v) in enumerate(edges):
        matrix[u,i] -= 1
        matrix[v,i] += 1
    # Construct a fundamental-cycle matrix with identity on the chords.
    reduced = Matrix(matrix[:-1].tolist())
    tree = reduced.rref()[1]
    chords = [i for i in range(len(edges)) if i not in tree]
    basis = np.zeros((len(edges),len(chords)),dtype=int)
    block = -reduced[:,tree].inv()*reduced[:,chords]
    assert all(x.q == 1 for x in block)
    basis[list(tree),:] = np.array(block.tolist(),dtype=int)
    basis[chords,:] = np.eye(len(chords),dtype=int)
    assert np.count_nonzero(matrix@basis) == 0
    minors = Counter()
    for rows in combinations(range(len(edges)),4):
        if len(chords) == 4:
            minors[int(Matrix(basis[list(rows)].tolist()).det())] += 1
        else:
            assert len(chords) < 4
            minors[0] += 1
    return dict(edges=edges, cycle_rank=len(chords), chords=chords,
                four_insertions={str(k):v for k,v in sorted(minors.items())},
                nonzero_pairs_of_placements=(sum(n for k,n in minors.items() if k))**2)


def sparse_cubic():
    data = json.loads(gzip.decompress((HERE/"representation26.json.gz").read_bytes()))
    rep = native26()
    weights = rep[:4].diagonal(axis1=1,axis2=2).T
    for k in range(4):
        assert np.array_equal(np.diag(weights[:,k]),rep[k])
    triples = list(combinations_with_replacement(range(26),3))
    unknown = [t for t in triples if not weights[list(t)].sum(axis=0).any()]
    lookup = {key:i for i,key in enumerate(unknown)}
    equations = set()
    for generator in rep[4:]:
        action = [list(zip(np.flatnonzero(col),col[col != 0])) for col in generator.T]
        for key in triples:
            row = Counter()
            for k,i in enumerate(key):
                for j, val in action[i]:
                    new = tuple(sorted((*key[:k],int(j),*key[k+1:])))
                    if new in lookup:
                        row[lookup[new]] += int(val)
            packed = tuple(sorted((i,x) for i,x in row.items() if x))
            if packed:
                equations.add(packed)
    equations = sorted(equations)
    mat = np.zeros((len(equations),len(unknown)),dtype=int)
    for i,row in enumerate(equations):
        for j,v in row:
            mat[i,j] = v
    kernel = Matrix(mat.tolist()).nullspace()
    assert len(kernel) == 1
    den = lcm(*(int(v.q) for v in kernel[0]))
    coefficients = [int(v*den) for v in kernel[0]]
    assert gcd(*coefficients) == 1
    cubic = np.zeros((26,)*3,dtype=np.int64)
    for key,x in zip(unknown,coefficients):
        for p in set(permutations(key)):
            cubic[p] = x
    for g in rep:
        # Full tensor Ward equation, including repeated slots and zero weights.
        value = np.einsum("ri,rjk->ijk",g,cubic)
        assert not np.any(value+value.transpose(1,2,0)+value.transpose(2,0,1))
    metric = Matrix.zeros(26)
    for i,j,val in data["form"]:
        metric[i,j] = Fraction(val)
    inverse = metric.inv()
    scale = lcm(*(int(val.q) for val in inverse))
    raw = np.array((scale*inverse).tolist(),dtype=np.int64)
    upper = np.einsum("ai,bj,ck,ijk->abc",raw,raw,raw,cubic,optimize=True)
    assert max(abs(int(x)) for x in cubic.flat)*26**3*max(abs(int(x)) for x in raw.flat)**3 < 2**62
    lower = {tuple(map(int,k)):int(cubic[tuple(k)]) for k in np.argwhere(cubic)}
    allowed = defaultdict(list)
    for i,j,k in np.argwhere(upper):
        allowed[int(i),int(j)].append((int(k),int(upper[i,j,k])))
    return weights,lower,allowed,dict(zero_weight_unknowns=len(unknown), Ward_rank=len(unknown)-1,
                                       full_Ward_entries=52*26**3, distinct_equations=len(equations),
                                       nonzero_lower=len(lower), inverse_metric_denominator=scale)


def graph_value(chords):
    weights,lower,allowed,record = sparse_cubic()
    @cache
    def determinant(keys, size):
        return int(Matrix(weights[list(keys),:size].tolist()).det())
    total = 0
    terms = Counter()
    for x, c in lower.items():
        for y, d in lower.items():
            for zs in product(*(allowed[x[i],y[i]] for i in range(3))):
                z = tuple(k for k,_ in zs)
                if z not in lower:
                    continue
                coeff = c*d*lower[z]*zs[0][1]*zs[1][1]*zs[2][1]
                total += coeff
                edge_labels = (*x,*y,*z)
                chosen = tuple(edge_labels[i] for i in chords)
                terms[chosen] += coeff
    assert total == sum(terms.values())
    pure = full = 0
    for key,coefficient in terms.items():
        if coefficient and len(set(key)) == 4:
            weight = determinant(key,4)**2
            if weight:
                pure += coefficient*weight
                full += coefficient*weight*determinant(key[:3],3)**2
    den = record["inverse_metric_denominator"]**9
    scores = [Fraction(x,den) for x in (total,pure,full)]
    return dict(cubic=record, raw_evaluations=list(map(str,scores)),
                # The ratio removes the choice of vector in the invariant cubic line.
                relative_to_undecorated=[str(x/scores[0]) for x in scores[1:]] if scores[0] else None,
                ratio_of_flag_to_pure=str(scores[2]/scores[1]) if scores[1] else None)


def calculate():
    tetrahedron = incidence(4,list(combinations(range(4),2)))
    bipartite = incidence(6,[(i,3+j) for i,j in product(range(3),repeat=2)])
    assert tetrahedron["nonzero_pairs_of_placements"] == 0
    assert bipartite["nonzero_pairs_of_placements"] == 81**2
    return dict(tetrahedron=tetrahedron, K33=bipartite,
                native_K33_Cartan=graph_value(bipartite["chords"]))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
