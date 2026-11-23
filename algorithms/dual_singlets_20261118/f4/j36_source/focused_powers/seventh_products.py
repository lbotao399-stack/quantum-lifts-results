"""Test the terminal six-colour tensors and the bilinear factors of the seventh source."""

from functools import cache
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import gzip
import json

import numpy as np
from flint import nmod_mat

from source_sizes import multiplicity, trace


HERE = Path(__file__).resolve().parent


@cache
def generators():
    path = HERE.parent / "global_boundaries" / "representation26.json.gz"
    entries = json.loads(gzip.decompress(path.read_bytes()))["matrices"]
    result = np.zeros((52, 26, 26), dtype=np.int64)
    for i, matrix in enumerate(entries):
        for j, k, value in matrix:
            result[i, j, k] = value
    diagonal = result[np.all(result == result * np.eye(26, dtype=int), axis=(1, 2))]
    assert len(diagonal) == 4
    return result, diagonal


def pairings(labels):
    if not labels:
        yield ()
        return
    i, *rest = labels
    for pos, j in enumerate(rest):
        for pairs in pairings(rest[:pos] + rest[pos + 1 :]):
            yield ((i, j), *pairs)


def six_colours(seed, prime, cartan=False):
    basis = generators()[int(cartan)]
    rng = np.random.Generator(np.random.PCG64(seed))
    coefficients = rng.integers(-4, 5, size=(6, len(basis)))
    colours = (coefficients @ basis.reshape(len(basis), -1)).reshape(6, 26, 26) % prime
    inverse_six = pow(6, -1, prime)

    def bracket(a, b):
        return (a @ b - b @ a) % prime

    def form(a, b):
        return (int((a * b.T).sum()) * inverse_six) % prime

    def symmetric_trace():
        # Fix the leading matrix by cyclicity. Symmetrize the remaining five.
        @cache
        def ordered_sum(mask):
            if not mask:
                return np.eye(26, dtype=np.int64)
            return sum(
                colours[j + 1] @ ordered_sum(mask ^ (1 << j))
                for j in range(5) if mask & (1 << j)
            ) % prime

        return int((colours[0] * ordered_sum(31).T).sum()) * pow(factorial(5), -1, prime) % prime

    products = []
    for pairs in pairings(list(range(6))):
        value = 1
        for i, j in pairs:
            value = value * form(colours[i], colours[j]) % prime
        products.append(value)
    products.append(symmetric_trace())
    if cartan:
        assert all(not bracket(colours[i], colours[j]).any() for i, j in combinations(range(6), 2))
        return products

    boundaries = []
    for i, j in combinations(range(6), 2):
        marked = [bracket(colours[i], colours[j]), *(colours[k] for k in range(6) if k not in (i, j))]
        for a, b in combinations(range(5), 2):
            c, d, e = (k for k in range(5) if k not in (a, b))
            boundaries.append(form(marked[a], marked[b]) * form(bracket(marked[c], marked[d]), marked[e]) % prime)
        for a, b, c in permutations((1, 2, 3)):
            nested = bracket(bracket(bracket(marked[0], marked[a]), marked[b]), marked[c])
            boundaries.append(form(nested, marked[4]))
    return boundaries + products


def pivots(matrix):
    reduced, rank = matrix.rref()
    return [next(j for j in range(reduced.ncols()) if reduced[i, j]) for i in range(rank)]


def minor(rows, columns, prime):
    return int(nmod_mat([[r[j] for j in columns] for r in rows], prime).det())


def remove_pair(shape):
    for i, j in combinations(range(len(shape)), 2):
        shortened = tuple(n - int(k in (i, j)) for k, n in enumerate(shape))
        if shortened == tuple(sorted(shortened, reverse=True)):
            yield tuple(n for n in shortened if n)


def first_target(shape):
    if len(shape) < 3 or shape[2] < 3:
        return None
    assert shape[2] == 3
    return (shape[0] + 1, shape[1] + 1, *shape[3:])


def run():
    dimension = trace((1,) * 6)
    assert dimension == 80
    minors = []
    for prime in (1009, 1013):
        samples = [six_colours(909737 + i, prime) for i in range(dimension)]
        cartan = [six_colours(316103 + i, prime, True) for i in range(16)]
        if prime == 1009:
            boundary_columns = pivots(nmod_mat([row[:240] for row in samples], prime))
            columns = boundary_columns + list(range(240, 256))
        assert len(boundary_columns) == 64
        full_det = minor(samples, columns, prime)
        cartan_det = int(nmod_mat(cartan, prime).det())
        assert full_det and cartan_det
        minors.append(dict(prime=prime, full_minor=full_det, cartan_minor=cartan_det))

    factors = []
    for shape in remove_pair((5, 5, 3, 3)):
        image = first_target(shape)
        factors.append(dict(source=shape, first_target=image, terminal_dimension=
                            multiplicity(image, sum(image), trace) if image else None))
    return dict(labelled_colour_invariants=dimension, bracket_insertions=240,
                independent_bracket_insertions=len(boundary_columns),
                metric_products=15, polarized_sextics=1, cartan_rank=16,
                boundary_columns=boundary_columns, ordered_minors=minors,
                bilinear_complements=factors)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
