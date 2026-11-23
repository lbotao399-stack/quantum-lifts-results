#!/usr/bin/env python3
"""Shrinking-dual updates on all actual SU3 active-profile Q0 rows.

This is an algorithm benchmark on a projected target, not a fortuitous census.
The same exact update applies when the streamed rows are full relation images.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np
from flint import nmod_mat

HERE = Path(__file__).resolve().parent


def matrix(a, prime):
    return nmod_mat(np.asarray(a, dtype=np.int64).tolist(), prime)


def kernel_columns(rows, prime):
    """Canonical columns spanning the exact annihilator of these rows."""
    rows = np.asarray(rows, dtype=np.int64)
    if rows.shape[1] == 0:
        return np.empty((0, 0), dtype=np.int64), 0
    reduced, rank = matrix(rows, prime).rref()
    pivots = []
    for i in range(rank):
        pivots.append(next(j for j in range(rows.shape[1]) if reduced[i, j]))
    free = [j for j in range(rows.shape[1]) if j not in pivots]
    dual = np.zeros((rows.shape[1], len(free)), dtype=np.int64)
    for k, j in enumerate(free):
        dual[j, k] = 1
        for i, pivot in enumerate(pivots):
            dual[pivot, k] = -int(reduced[i, j]) % prime
    assert not np.any(rows @ dual % prime)
    return dual, rank


def update_duals(duals, new_rows, prime, identity=False):
    """Ann(A+B) = Ann(A) ker(B Ann(A)); no source nullspace is formed."""
    residual = new_rows if identity else new_rows @ duals % prime
    residual_duals, added_rank = kernel_columns(residual, prime)
    new_duals = residual_duals if identity else duals @ residual_duals % prime
    return new_duals, added_rank, residual.shape


def one(source, prime):
    suffix = '' if prime == 1000033 else '_p1000003'
    directory = source / ('projected_Q0' + suffix)
    cert_path = directory / 'projected_Q0_certificate.json'
    data_path = directory / ('projected_Q0_basis_coordinates_p%d.npy' % prime)
    certificate = json.loads(cert_path.read_text())
    rows = np.load(data_path).astype(np.int64)
    ids = np.asarray(certificate['source_row_ids'])
    n = rows.shape[1]
    blocks = [rows[ids[:, 0] == i] for i in (38, 39, 41, 53)]

    def incremental():
        duals = None
        profile_data = []
        total_rank = 0
        for i, (profile, block) in enumerate(zip((38, 39, 41, 53), blocks)):
            before = n if i == 0 else duals.shape[1]
            duals, added, shape = update_duals(duals, block, prime, identity=i == 0)
            total_rank += added
            profile_data.append({'source_profile': profile, 'source_rows': len(block),
                'dual_dimension_before': before, 'added_rank': added,
                'dual_dimension_after': duals.shape[1],
                'residual_matrix_shape': list(shape)})
        return duals, total_rank, profile_data

    trials = []
    # Both methods read the same already generated dense rows; no generation
    # cost is included. Alternate order to reduce systematic timing bias.
    for trial in range(4):
        sequence = ('full', 'incremental') if trial % 2 == 0 else ('incremental', 'full')
        times = {}
        for method in sequence:
            start = time.perf_counter_ns()
            if method == 'full':
                full_duals, full_rank = kernel_columns(rows, prime)
            else:
                duals, rank, profile_data = incremental()
            times[method] = time.perf_counter_ns() - start
        assert full_rank == rank == 430
        assert duals.shape == full_duals.shape == (434, 4)
        assert matrix(np.column_stack((duals, full_duals)), prime).rank() == 4
        trials.append(times)
    assert not np.any(rows @ duals % prime)
    rhs_certificate = json.loads((source / ('active_frame' + suffix) / 'active_second_prime_certificate.json').read_text())
    rhs = np.asarray(rhs_certificate['RHS_basis_coefficients'], dtype=np.int64)
    assert not np.any(rhs @ duals % prime)
    return {'prime': prime, 'source_rows': len(rows), 'target_dimension': n,
        'rank': rank, 'source_kernel_dimension_not_constructed': len(rows) - rank,
        'final_dual_dimension': duals.shape[1], 'profiles': profile_data,
        'full_elimination_matrix_entries': int(rows.size),
        'incremental_elimination_matrix_entries': sum(x['residual_matrix_shape'][0] * x['residual_matrix_shape'][1] for x in profile_data),
        'all_source_dual_pairings_checked': int(rows.shape[0] * duals.shape[1]),
        'RHS_dual_pairings': (rhs @ duals % prime).tolist(),
        'dual_columns': duals.tolist(), 'timings_ns': trials,
        'input_sha256': {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [cert_path, data_path]}}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=HERE.parents[1] / 'quantum_lifts_frontier_work/su3_second')
    ap.add_argument('--output', type=Path, default=HERE / 'streaming_profile_duals.json')
    args = ap.parse_args()
    result = {'passed': True, 'runs': [one(args.input, p) for p in (1000003, 1000033)],
        'scope': 'Exact modular algorithm control on the complete projected target31 image only. Not a full coupled second correction, cohomology calculation or fortuitous census. Timings include arithmetic on existing arrays, not native differential or trace generation. The update is ordinary exact linear algebra; the contribution is its checked application and exported complete duals.'}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'passed': True, 'runs': [{k: v for k, v in x.items() if k not in ('dual_columns', 'input_sha256')} for x in result['runs']]}, indent=2))


if __name__ == '__main__':
    main()
