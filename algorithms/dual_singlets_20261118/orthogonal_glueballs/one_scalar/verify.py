"""Check the second-page dual on the complete finite-colour one-scalar relations."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

from flint import fmpq_mat, nmod_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from orthogonal_one_scalar import dimension, dual, relation, specification
from terminal_commutators import words


def pivots(matrix):
    reduced, rank = matrix.rref()
    return [next(j for j in range(reduced.ncols()) if reduced[i, j]) for i in range(rank)]


def residues(values, prime):
    return [v.numerator * pow(v.denominator, -1, prime) % prime for v in map(Fraction, values)]


def check_row(task):
    row, n, k, probes, plane_indices, sites, scores, prime = task
    assert row['spec'] == specification(row['index'], n, k)
    pfaffian = relation(k, row['spec'])
    assert pfaffian.pf_degree > n
    ordered = [pfaffian.coefficient(w) for w in probes]
    for j, site in zip(plane_indices, sites, strict=True):
        assert pfaffian.cartan(site) == Fraction(row['cartan'][j])
    assert sum(c * pfaffian.coefficient(w) for w, c in scores) == Fraction(row['second_score'])
    pfaffian.clear()
    return residues(ordered, prime)


def verify(n, k, rebuild, workers):
    data = json.loads(gzip.decompress((HERE / f'one-scalar-relations-{n}-{k}.json.gz').read_bytes()))
    assert (data['colours'], data['k']) == (n, k)
    rows = data['rows']
    lower = [dimension(n, k, r, full_weight=True) for r in (0, 1)]
    count = dimension(n, k, 2)
    assert lower == data['lower'] == [0, 0]
    assert count == data['relation_dimension'] == data['ordered_rank']
    assert all(relation(k, r['spec']).pf_degree > n for r in rows)

    prime = 1000003
    ordered = nmod_mat([residues(row['ordered'], prime) for row in rows], prime)
    spanning = pivots(ordered.transpose())
    cols = pivots(ordered)
    assert len(spanning) == len(cols) == count
    minor = [[int(ordered[i, j]) for j in cols] for i in spanning]
    det = int(nmod_mat(minor, prime).det())
    assert det

    planes = fmpq_mat([row['cartan'] for row in rows])
    values = fmpq_mat([[row['second_score']] for row in rows])
    plane_indices = pivots(planes)
    plane_rows = pivots(planes.transpose())
    square = fmpq_mat([[planes[i, j] for j in plane_indices] for i in plane_rows])
    weights = square.solve(fmpq_mat([[values[i, 0]] for i in plane_rows]))
    selected = fmpq_mat([[planes[i, j] for j in plane_indices] for i in range(len(rows))])
    assert selected * weights == values
    assert len(plane_indices) == data['cartan_rank'] == data['augmented_rank']

    adjoint = dual(n, k)
    assert [adjoint[grade]['checked_columns'] for grade in ('first', 'second')] == [55, 2178]
    if rebuild:
        probes = [tuple(map(tuple, data['probes'][j])) for j in cols]
        all_sites = words((k-2, k-3, 3, 2, 2))
        sites = [all_sites[j] for j in plane_indices]
        scores = [(tuple(map(tuple, w)), Fraction(c)) for w, c in adjoint['scores']]
        tasks = ((rows[i], n, k, probes, plane_indices, sites, scores, prime) for i in spanning)
        with ProcessPoolExecutor(workers) as pool:
            for expected, actual in zip(minor, pool.map(check_row, tasks), strict=True):
                assert actual == expected

    return dict(colours=n, k=k, lower_relations=lower, second_relations=count,
                determinant_mod_1000003=det, plane_indices=plane_indices,
                weights=[str(v[0]) for v in weights.tolist()], nonzero_residuals=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild', action='store_true', help='recompute every basis Pfaffian')
    parser.add_argument('--workers', type=int, default=3)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error('enable Python assertions and give a positive worker count')
    for n, k in [(11, 6), (13, 7)]:
        print(json.dumps(verify(n, k, args.rebuild, args.workers)), flush=True)
