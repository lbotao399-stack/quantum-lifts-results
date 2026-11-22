"""Split G2 over the rationals in its seven-dimensional fundamental module."""

from functools import cache
from math import gcd

from flint import fmpq_mat, fmpz_mat


@cache
def structure():
    from bmn_primary_index import G2_ROOTS

    weights = ((1, 0), (-1, 1), (2, -1), (0, 0), (-2, 1), (1, -1), (-1, 0))
    simple = ((2, -1), (-3, 2))
    basis = []
    for root in simple:
        raising, lowering = fmpz_mat(7, 7), fmpz_mat(7, 7)
        for i, (x, y) in enumerate(weights):
            if (x + root[0], y + root[1]) in weights:
                continue
            length = 0
            while (x - (length + 1) * root[0], y - (length + 1) * root[1]) in weights:
                length += 1
            previous = i
            for k in range(1, length + 1):
                j = weights.index((x - k * root[0], y - k * root[1]))
                raising[previous, j] = length - k + 1
                lowering[j, previous] = k
                previous = j
        basis.append((root, raising, lowering))

    def bracket(a, b):
        return a * b - b * a

    cartan = [bracket(e, f) for _, e, f in basis]
    for j, h in enumerate(cartan):
        assert h == fmpz_mat([[w[j] * int(i == k) for k in range(7)] for i, w in enumerate(weights)])
    root_vectors = {}
    for sign in (1, -1):
        simple_vectors = [(tuple(sign * i for i in root), e if sign == 1 else f) for root, e, f in basis]
        pending = list(simple_vectors)
        root_vectors.update(simple_vectors)
        while pending:
            weight, matrix = pending.pop()
            for root, generator in simple_vectors:
                total = tuple(a + b for a, b in zip(weight, root))
                if total not in G2_ROOTS or total in root_vectors:
                    continue
                result = bracket(matrix, generator)
                divisor = gcd(*map(int, result.entries()))
                assert divisor
                result = fmpz_mat(7, 7, [int(v) // divisor for v in result.entries()])
                root_vectors[total] = result
                pending.append((total, result))
    assert set(root_vectors) == set(G2_ROOTS)
    roots = [(0, 0), (0, 0), *G2_ROOTS]
    matrices = cartan + [root_vectors[r] for r in G2_ROOTS]
    coordinates = fmpq_mat([m.entries() for m in matrices]).transpose()
    left = (coordinates.transpose() * coordinates).inv() * coordinates.transpose()
    constants = {}
    for i, a in enumerate(matrices):
        for j, b in enumerate(matrices):
            column = fmpq_mat(49, 1, bracket(a, b).entries())
            coefficients = left * column
            assert coordinates * coefficients == column
            terms = []
            for k, value in enumerate(coefficients.entries()):
                if value:
                    assert value.denominator == 1
                    assert roots[k] == tuple(x + y for x, y in zip(roots[i], roots[j]))
                    terms.append((k, int(value)))
            constants[i, j] = terms
    return roots, constants
