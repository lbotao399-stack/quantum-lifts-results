"""Reuse universal trace-product rows on an explicitly certified G2 torus.

The seven eigenvalue weights are x,y,-x-y,-x,-y,x+y,0. They preserve
a stable three-form. Its GL7 stabilizer has dimension 14 by rational
linear algebra. Dropping the final zero eigenvalue changes no trace.
"""
from collections import defaultdict
from itertools import combinations, permutations
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import signal
import time

import sp_direct_hw_cartan_products_v1 as direct
import numpy as np
from flint import fmpq_mat, nmod_mat


def torus_certificate():
    weights = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1), (0, 0)]
    triples = [(0, 1, 2), (3, 4, 5), (0, 3, 6), (1, 4, 6), (2, 5, 6)]
    phi = np.zeros((7, 7, 7), dtype=np.int64)
    for triple in triples:
        assert all(sum(weights[index][component] for index in triple) == 0 for component in (0, 1))
        for order in permutations(range(3)):
            sign = (-1) ** sum(order[i] > order[j] for i in range(3) for j in range(i + 1, 3))
            phi[tuple(triple[index] for index in order)] = sign
    equations = []
    for i, j, k in combinations(range(7), 3):
        row = []
        for r in range(7):
            for s in range(7):
                row.append(int((phi[r, j, k] if s == i else 0)
                               + (phi[i, r, k] if s == j else 0)
                               + (phi[i, j, r] if s == k else 0)))
        equations.append(row)
    rank = int(fmpq_mat(equations).rank())
    assert rank == 35
    for component in (0, 1):
        diagonal = np.diag([weight[component] for weight in weights])
        assert not np.asarray(equations).dot(diagonal.reshape(49)).any()
    return dict(status='rational_stable_three_form_torus_certificate',
                eigenvalue_weights=weights, nonzero_oriented_triples=triples,
                GL7_dimension=49, infinitesimal_orbit_rank=rank, stabilizer_dimension=49 - rank,
                diagonal_stabilizer_checks=2,
                interpretation='A three-form with open GL7 orbit has connected complex stabilizer G2. These two independent diagonal matrices lie in its Lie algebra. Their eigenvalues span the rank-two torus; the trace in seven dimensions equals the six nonzero eigenvalue trace.',
                equation_matrix_sha256=hashlib.sha256(np.asarray(equations, dtype=np.int64).tobytes()).hexdigest())


def restrict_frame(frame):
    direct.diagonalize(frame)
    for letter, (odd, values) in frame.L.items():
        values[..., 2, 2] = (-values[..., 0, 0] - values[..., 1, 1]) % frame.p
        for index in range(3):
            values[..., index + 3, index + 3] = -values[..., index, index] % frame.p
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=int, default=240)
    args = parser.parse_args()
    start = time.monotonic()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('G2 Cartan reuse cap')))
    signal.alarm(args.seconds)
    source = json.loads((args.input / 'result.json').read_text())
    assert source['status'] == 'direct_highest_weight_Cartan_rank_lower_bound'
    assert direct.sha(args.input / 'family.pkl') == source['family_sha256']
    family = pickle.loads((args.input / 'family.pkl').read_bytes())
    prime = source['prime']
    direct.sp_census.configure(3, prime)
    certificate = torus_certificate()
    # No finite-rank relations enter the product family or its raising kernel.
    for generator in family['generators']:
        assert not direct.sp_census.base.action_expr(generator['expr'])
    kernel = family['kernel']
    assert int(max(np.diff(kernel.indptr), default=0)) * (prime - 1) ** 2 < 2 ** 63
    blocks, comparisons = [], 0
    for metadata in source['frames']:
        degree, points, seed, k = [metadata[key] for key in ('degree', 'points', 'seed', 'k')]
        frame = restrict_frame(direct.sp_census.base.Frame(degree, points, seed, prime, k=k))
        necklaces = {traces[0] for generator in family['generators'] for paths, traces in generator['expr']}
        probes = []
        for d in range(min(k, 2) + 1):
            candidates = sorted(trace for trace in necklaces if sum(3 <= letter <= 5 for letter in trace) == d)
            if candidates:
                probes.extend(candidates[index] for index in sorted({0, len(candidates) // 2, len(candidates) - 1}))
        actual = direct.diagonal_necklaces(frame, probes)
        expected = frame.necklace_values(probes)
        for trace in probes:
            assert actual[trace][1] == expected[trace][1]
            assert np.array_equal(actual[trace][0], expected[trace][0])
            comparisons += actual[trace][0].size
        frame.necklace_values = lambda words, frame=frame: direct.diagonal_necklaces(frame, words)
        raw = direct.product_values(family['generators'], family['products'], frame, degree)
        blocks.append(np.asarray(kernel @ raw) % prime)
    values = np.concatenate(blocks, axis=1)
    selected = direct.independent_rows(values, prime) if values.size else []
    if selected:
        reduced, rank = nmod_mat(values[selected].tolist(), prime).rref()
        assert int(rank) == len(selected)
        pivots = [next(index for index, value in enumerate(row) if value) for row in reduced.tolist()[:int(rank)]]
        determinant = int(nmod_mat(values[np.ix_(selected, pivots)].tolist(), prime).det())
        assert determinant
    else:
        pivots, determinant = [], 1
    args.output.mkdir(parents=True, exist_ok=True)
    np.save(args.output / 'cartan_values.npy', values, allow_pickle=False)
    result = dict(status='direct_highest_weight_G2_Cartan_rank_lower_bound', group='G2',
                  sector=source['sector'], prime=prime, rank_lower_bound=len(selected),
                  universal_source_family_sha256=source['family_sha256'],
                  source_result_sha256=direct.sha(args.input / 'result.json'),
                  frames=source['frames'], torus_certificate=certificate,
                  independent_matrix_trace_coordinates=int(comparisons),
                  minor_rows=list(map(int, selected)), minor_columns=pivots, minor_determinant=determinant,
                  values_sha256=direct.sha(args.output / 'cartan_values.npy'),
                  driver_sha256=direct.sha(__file__), seconds=time.monotonic() - start,
                  scope='Direct highest-weight lower bound for the complexified G2 gauge algebra, hence over Q. Universal formal trace rows are reused without any Sp3 physical-quotient relation. No complete Cartan image or quantum obstruction is asserted.')
    direct.save(args.output / 'result.json', result)
    print(json.dumps({key: value for key, value in result.items() if key not in ('minor_rows', 'minor_columns', 'frames')}), flush=True)


if __name__ == '__main__':
    main()
