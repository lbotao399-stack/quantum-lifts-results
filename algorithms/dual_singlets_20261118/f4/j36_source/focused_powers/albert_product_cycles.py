"""Young projection of the exceptional nine-scalar times two invariant bilinears."""
import argparse
import json
from fractions import Fraction
from math import lcm
from concurrent.futures import ProcessPoolExecutor
from functools import cache, partial
from itertools import combinations, combinations_with_replacement, permutations, product
from multiprocessing import get_context
from pathlib import Path

import numpy as np

from albert_boundaries import PAIRS, action, apply, balanced, colour_tensors, samples

AXES = tuple(range(4))
OMITTED = tuple(combinations_with_replacement(AXES, 3))
GRAMS = tuple(combinations_with_replacement(AXES, 2))
WEIGHT = (4, 4, 4, 2)


def sign(labels):
    return (-1) ** sum(a > b for i, a in enumerate(labels) for b in labels[i+1:])


def weight(column):
    omitted, gram, pair, last = column
    return tuple(3 - omitted.count(i) + gram.count(i) + pair.count(i) + int(i == last) for i in AXES)


def raised(column, i):
    omitted, gram, pair, last = column
    j = i + 1
    if i in omitted:
        v = list(omitted)
        v[v.index(i)] = j
        yield -omitted.count(i), (tuple(sorted(v)), gram, pair, last)
    if j in gram:
        v = list(gram)
        v[v.index(j)] = i
        yield gram.count(j), (omitted, tuple(sorted(v)), pair, last)
    if j in pair and i not in pair:
        v = tuple(i if n == j else n for n in pair)
        yield sign(v), (omitted, gram, tuple(sorted(v)), last)
    if last == j:
        yield 1, (omitted, gram, pair, i)


@cache
def highest():
    columns = [t for t in product(OMITTED, GRAMS, PAIRS, AXES) if weight(t) == WEIGHT]
    relations = {}
    def add(key, col, value):
        row = relations.setdefault(key, [0]*len(columns))
        row[col] += value
    for col, t in enumerate(columns):
        for i in range(3):
            for s,v in raised(t, i):
                add((i,v), col, s)
        omitted, gram, pair, k = t
        if k not in pair:
            ordered = (*pair,k)
            add(('split',omitted,gram,tuple(sorted(ordered))), col, sign(ordered))
    matrix = [[Fraction(v) for v in row] for row in relations.values()]
    leading = []
    for j in range(len(columns)):
        r = len(leading)
        pivot = next((i for i in range(r, len(matrix)) if matrix[i][j]), None)
        if pivot is None:
            continue
        matrix[r], matrix[pivot] = matrix[pivot], matrix[r]
        row = [v / matrix[r][j] for v in matrix[r]]
        matrix[r] = row
        for i, old in enumerate(matrix):
            if i != r and old[j]:
                a = old[j]
                matrix[i] = [x - a*y for x,y in zip(old,row)]
        leading.append(j)
    vectors = []
    for k in range(len(columns)):
        if k in leading:
            continue
        vector = [Fraction(0)] * len(columns)
        vector[k] = Fraction(1)
        for i,j in enumerate(leading):
            vector[j] = -matrix[i][k]
        denominator = lcm(*(v.denominator for v in vector))
        vectors.append([int(v * denominator) for v in vector])
    coeff = np.array(vectors, dtype=np.int64)
    assert coeff.shape == (1,96)
    assert not np.any(np.asarray(list(relations.values()),dtype=np.int64) @ coeff.T)
    return columns, coeff


def alternating(corner, actions, omitted, p):
    indices = tuple(i for i in AXES if i != omitted)
    out = np.zeros_like(corner)
    for order in permutations(indices):
        term = corner
        for axis, i in enumerate(order):
            term = apply(actions[i], term, (axis,))
            if axis == 0:
                term = balanced(term, p)
        out += sign(order)*term
    return (-1)**omitted * balanced(out, p)


def scalar_polynomials(x, p, active):
    lower, upper = colour_tensors()[:2]
    up = [action(t)[0] for t in x]
    down = [action(t.T)[0] for t in x]
    ids = sorted({i for t in active for i in t})
    right = {i: alternating(upper, up, i, p) for i in ids}
    left = np.stack([alternating(lower, down, i, p).ravel() for i in ids])
    triples = {}
    for j,k in sorted({(j,k) for t in active for _,j,k in permutations(t)}):
        mid = alternating(right[k], up, j, p)
        assert 26**4*(p//2)**2 < 2**53
        scores = left @ mid.ravel()
        assert np.array_equal(scores, scores.astype(np.int64))
        triples.update({(i,j,k): int(v)%p for i,v in zip(ids,scores)})
    return {t: sum(triples[o] for o in permutations(t)) % p for t in active}


def point(n, p, columns, coeff, closed=False):
    x, _, y = samples(n, p, closed)
    scalars = scalar_polynomials(x, p, {t[0] for t in columns})
    gram = x.reshape(4,-1) @ x.transpose(0,2,1).reshape(4,-1).T % p
    marked = y.reshape(6,-1) @ x.transpose(0,2,1).reshape(4,-1).T % p
    v = [scalars[t] * gram[g] * marked[PAIRS.index(pair),k] % p for t,g,pair,k in columns]
    return np.asarray((coeff @ v)%p, dtype=np.uint16)


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--prime', type=int, default=1009)
    a.add_argument('--workers', type=int, default=1)
    a.add_argument('--points', type=Path)
    a.add_argument('--output', type=Path)
    args = a.parse_args()
    cols, coeff = highest()
    keep = np.flatnonzero(np.any(coeff,axis=0))
    cols, coeff = [cols[i] for i in keep], coeff[:,keep]
    print('formal weights', len(highest()[0]), 'highest closed',len(coeff),'active',len(cols),flush=True)
    f=partial(point, p=args.prime, columns=cols, coeff=coeff)
    vals=f(32)
    assert vals.any() and not f(32,closed=True).any()
    print('free and split',vals.tolist(),[0]*len(vals),flush=True)
    if args.points:
        ids=json.loads(args.points.read_bytes())
        with ProcessPoolExecutor(args.workers, mp_context=get_context('fork')) as pool:
            out=np.stack(list(pool.map(f,ids,chunksize=2)))
        np.savez_compressed(args.output,points=ids,values=out)


if __name__ == '__main__':
    main()
