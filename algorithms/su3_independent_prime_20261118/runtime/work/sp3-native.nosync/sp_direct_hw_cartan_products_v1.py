"""Rational lower bounds from direct highest-weight Cartan product ranks.

No subtraction of uncertified weight-space ranks occurs. Each modular
highest-weight product row admits the usual rational extremal projection.
Cartan restriction kills every Q0 boundary, so a nonzero minor supplies
independent rational classical classes.
"""
from collections import defaultdict
from pathlib import Path
from math import comb
import argparse
import hashlib
import json
import pickle
import signal
import time

import sp_census
import numpy as np
from flint import nmod_mat
from fortuity.bmn.words import graviton_generators
from fortuity.hwcoh.evalg import table
from protected import graph_expr, expr_sec, monomials
from hw_products import Raising
from sp_graviton_singlets import product_values
from sp_sparse_formal_kernel_v1 import split_kernel_sparse
from sp_projected_basis import independent_rows

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def diagonalize(frame):
    for letter, (odd, values) in list(frame.L.items()):
        diagonal = np.zeros_like(values)
        for index in range(6):
            diagonal[..., index, index] = values[..., index, index]
        assert np.array_equal(diagonal[..., :3, :3], -diagonal[..., 3:, 3:] % frame.p)
        frame.L[letter] = (odd, diagonal)
    return frame


def diagonal_necklaces(frame, necklaces):
    result = {}
    for trace in necklaces:
        degree = sum(3 <= letter <= 5 for letter in trace)
        if degree > frame.k:
            result[trace] = np.zeros((frame.P, 0), dtype=np.int64), degree
            continue
        current = np.ones((frame.P, 1, 6), dtype=np.int64)
        used = 0
        for letter in trace:
            odd, matrix = frame.L[letter]
            diagonal = np.diagonal(matrix, axis1=-2, axis2=-1)
            if not odd:
                current = current * diagonal % frame.p
            else:
                left, right, signs, starts = table(frame.k, used, 1)
                current = np.add.reduceat((current[:, left] * diagonal[:, right] % frame.p)
                                          * signs[None, :, None], starts, axis=1) % frame.p
                used += 1
        assert used == degree
        result[trace] = current.sum(axis=-1) % frame.p, degree
    return result


def generators(prime):
    values = []
    for name, expression in graviton_generators(6, 1):
        expression = graph_expr(expression)
        if not expression:
            continue
        assert not sp_census.base.action_expr(expression)
        lengths = {sum(map(len, traces)) for paths, traces in expression}
        assert len(lengths) == 1 and next(iter(lengths)) % 2 == 0
        odd = {sum(sp_census.base.count(word)[3:6]) % 2 for word in expression}
        assert len(odd) == 1
        values.append(dict(name=name, expr=expression, sector=expr_sec(expression),
                           odd=odd.pop(), reflection=0, family='trace'))
    return values


def build(sector, prime, points, output):
    start = time.monotonic()
    sp_census.configure(3, prime)
    family = generators(prime)
    products = monomials(family, sector, 0)
    raising = Raising(family, prime)
    images = [{(root, word): coefficient for root in (1, 2)
               for word, coefficient in raising.product(recipe, root).items()} for recipe in products]
    kernel, rank = split_kernel_sparse(images, prime)
    assert kernel.shape == (len(products) - rank, len(products))
    for row in range(kernel.shape[0]):
        residual = defaultdict(int)
        for offset in range(kernel.indptr[row], kernel.indptr[row + 1]):
            for word, coefficient in images[int(kernel.indices[offset])].items():
                residual[word] += int(kernel.data[offset]) * coefficient
        assert not any(coefficient % prime for coefficient in residual.values())
    output.mkdir(parents=True, exist_ok=True)
    (output / 'family.pkl').write_bytes(pickle.dumps(dict(generators=family, products=products, kernel=kernel), protocol=5))
    print(json.dumps(dict(stage='product_highest_weights', sector=sector, weight_products=len(products),
                          highest_weight_rows=kernel.shape[0], seconds=time.monotonic() - start)), flush=True)
    assert int(max(np.diff(kernel.indptr), default=0)) * (prime - 1) ** 2 < 2 ** 63
    blocks, frames, comparisons = [], [], 0
    degrees = range(sector[0] % 2, min(sector[0], 9) + 1, 2)
    for degree in degrees:
        k = min(9, degree + 2) if degree else 0
        seed = 2026111761 + 7919 * degree
        frame = diagonalize(sp_census.base.Frame(degree, points, seed, prime, k=k))
        necklaces = {traces[0] for generator in family for paths, traces in generator['expr']}
        # Compare the optimized scalar exterior products with the independent
        # dense-matrix trace evaluator at the same points before using them.
        probes = []
        for probe_degree in range(min(k, 2) + 1):
            candidates = sorted(trace for trace in necklaces
                                if sum(3 <= letter <= 5 for letter in trace) == probe_degree)
            if candidates:
                probes.extend(candidates[index] for index in sorted({0, len(candidates) // 2, len(candidates) - 1}))
        expected = frame.necklace_values(probes)
        actual = diagonal_necklaces(frame, probes)
        for trace in probes:
            assert expected[trace][1] == actual[trace][1]
            assert np.array_equal(expected[trace][0], actual[trace][0])
            comparisons += expected[trace][0].size
        frame.necklace_values = lambda words, frame=frame: diagonal_necklaces(frame, words)
        raw = product_values(family, products, frame, degree)
        block = np.asarray(kernel @ raw) % prime
        blocks.append(block)
        frames.append(dict(degree=degree, k=k, points=points, seed=seed, columns=frame.ncols))
    values = np.concatenate(blocks, axis=1) if blocks else np.zeros((kernel.shape[0], 0), dtype=np.int64)
    np.save(output / 'cartan_values.npy', values, allow_pickle=False)
    selected = independent_rows(values, prime) if values.size else []
    if selected:
        reduced, rank = nmod_mat(values[selected].tolist(), prime).rref()
        assert int(rank) == len(selected)
        pivots = [next(index for index, value in enumerate(row) if value) for row in reduced.tolist()[:int(rank)]]
        determinant = int(nmod_mat(values[np.ix_(selected, pivots)].tolist(), prime).det())
        assert determinant
    else:
        pivots, determinant = [], 1
    result = dict(status='direct_highest_weight_Cartan_rank_lower_bound', group='Sp(3)', sector=sector,
                  prime=prime, rank_lower_bound=len(selected), weight_products=len(products),
                  highest_weight_product_rows=kernel.shape[0], exact_closed_generators=len(family),
                  frames=frames, independent_matrix_trace_coordinates=int(comparisons),
                  minor_rows=list(map(int, selected)), minor_columns=pivots, minor_determinant=determinant,
                  family_sha256=sha(output / 'family.pkl'), values_sha256=sha(output / 'cartan_values.npy'),
                  driver_sha256=sha(__file__), seconds=time.monotonic() - start,
                  rational_bridge='Integer products are exactly Q0 closed. Apply the pinned rational extremal projector to integer lifts of modular product-kernel rows. It preserves closure and their mod-p reductions; its elementary denominators are below the prime. Cartan restriction kills Q0 boundaries.',
                  scope='A lower bound on independent rational closed products, obtained directly in highest weight. No signed difference of sampled ranks, cohomology completeness, or quantum lift is inferred.')
    save(output / 'result.json', result)
    print(json.dumps({key: value for key, value in result.items() if key not in ['minor_rows', 'minor_columns', 'frames']}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sector', required=True)
    parser.add_argument('--prime', type=int, default=1000003)
    parser.add_argument('--points', type=int, default=8)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seconds', type=int, default=600)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('direct Cartan product cap')))
    signal.alarm(args.seconds)
    build(tuple(map(int, args.sector.split(','))), args.prime, args.points, args.output)
