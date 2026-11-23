"""Check the four identities and the full rank of the 43 level-42 exchange recipes."""

import argparse
from fractions import Fraction
import gzip
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, eye

from seventh_crossings import PAIRS, point
from albert_boundaries import native26
from integer_section import FIRST_PAGE, HERE as SECTION, tensors


def native_identities():
    dimension = 26
    source = json.loads(gzip.decompress((FIRST_PAGE / 'representation26.json.gz').read_bytes()))
    form = np.zeros((dimension, dimension), dtype=np.int64)
    for i, j, value in source['form']:
        entry = Fraction(value)
        assert entry.denominator == 1
        form[i, j] = entry.numerator
    raised = 3 * form
    raised[:2, :2] = [[2, 1], [1, 2]]
    assert np.array_equal(raised @ form, 3 * np.eye(dimension, dtype=np.int64))
    cubic = np.zeros((dimension,) * 3, dtype=np.int64)
    for i, j, k, value in json.loads((SECTION / 'cubic.json').read_bytes())['cubic']:
        cubic[i, j, k] = value
    assert np.array_equal(cubic, cubic.transpose(1, 2, 0))
    assert np.array_equal(cubic, cubic.transpose(0, 2, 1))
    generators = native26().astype(np.int64)
    assert np.array_equal(form @ generators, -generators.transpose(0, 2, 1) @ form)
    products = generators[:, None] @ generators[None, :]
    jordan = raised @ cubic
    assert not np.any(products.reshape(52**2, -1) @ jordan.transpose(0, 2, 1).reshape(dimension, -1).T)
    _, upper, exchange, scales = tensors()
    assert scales == (28, 252, 84)
    raised_cubic = (raised @ cubic) @ raised
    generated = np.einsum('ab,bik,ajl->ijkl', raised, raised_cubic, raised_cubic, optimize=True)
    position = tuple(np.argwhere(upper)[0])
    assert np.array_equal(int(upper[position]) * generated, int(generated[position]) * upper)
    exchange = csr_matrix(exchange, dtype=np.int64)
    identity = eye(dimension**2, format='csr', dtype=np.int64)
    swap = identity[np.arange(dimension**2).reshape(dimension, dimension).T.ravel()]
    polynomial = (exchange - 6 * identity) @ (exchange + 36 * identity) @ (identity + swap)
    contraction = csr_matrix(np.outer(raised.ravel(), form.ravel()))
    assert (polynomial - 240 * contraction).nnz == 0
    assert (swap @ exchange - exchange @ swap).nnz == 0
    return {'native_generators': len(generators), 'ordered_jordan_trace_tests': 26*52**2,
            'two_slot_matrix_order': dimension**2, 'symmetric_trace_coefficient': 240}


def relation_columns():
    relations = []
    for i,j in combinations(range(len(PAIRS)), 2):
        if set(PAIRS[i]).isdisjoint(PAIRS[j]):
            row = [0] * 43
            row[7 + 6*i+j], row[7 + 6*j+i] = -1, 1
            relations.append(row)
    i = PAIRS.index((1, 3))
    row = [0] * 43
    row[0], row[1+i], row[7+6*i+i] = -216, 30, 1
    relations.append(row)
    return np.asarray(relations, dtype=np.int64)


def determinant(matrix, prime):
    rows = [[int(x) % prime for x in row] for row in matrix]
    result = 1
    for j in range(len(rows)):
        pivot = next(k for k in range(j, len(rows)) if rows[k][j])
        if pivot != j:
            rows[pivot], rows[j] = rows[j], rows[pivot]
            result = -result
        value = rows[j][j]
        result = result * value % prime
        inverse = pow(value, -1, prime)
        for k in range(j+1, len(rows)):
            m = rows[k][j] * inverse % prime
            for h in range(j, len(rows)):
                rows[k][h] = (rows[k][h] - m * rows[j][h]) % prime
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--regenerate', action='store_true', help='regenerate all native point evaluations')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    report = native_identities()
    relations = relation_columns()
    assert len(relations) == 4
    checks = []
    for prime in (1009, 1013):
        cert = json.loads((here / f'seventh-crossings-{prime}.json').read_bytes())
        with np.load(here / f'seventh-rows-{prime}.npz', allow_pickle=False) as data:
            rows, free = data['closed'].astype(np.int64), data['free'].astype(np.int64)
            assert data['point_ids'].tolist() == cert['points']
        if args.regenerate:
            assert np.array_equal(rows, np.stack([point(i, prime) for i in cert['points']]))
            assert np.array_equal(free, np.stack([point(i, prime, True).reshape(2,43) for i in (3,27)]))
        assert rows.shape == (39,43) and free.shape == (2,2,43)
        assert not np.any(rows @ relations.T % prime) and not np.any(free @ relations.T % prime)
        value = determinant(rows[:,cert['columns']], prime)
        assert value == cert['determinant']
        checks.append(dict(prime=prime, ordered_minor=value))
    report.update(exchange_recipes=43, left_and_right_recipes=86, independent_cochain_identities=4,
                  outgoing_rank=39, formal_cycle_dimension=47,
                  quotient_by_classical_boundaries_dimension=0, fixed_minors=checks)
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
