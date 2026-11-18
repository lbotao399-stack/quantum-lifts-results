"""Independent-node continuation of one source fermion degree.

Reuses the saved rational projection recipes and original evaluation frames.
Individual degree ranks must be combined before drawing cohomology conclusions.
"""
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
import argparse
import hashlib
import json
import multiprocessing
import os
import pickle
import time

import g2_bmn_complex_v1 as base
import numpy as np
from flint import nmod_mat
from extend_rows_binary64_exact_v1 import extend_rows
from sp_projected_basis import independent_rows


def reduce(values, prime):
    if not len(values):
        return values.copy(), []
    reduced, rank = nmod_mat(values.tolist(), prime).rref()
    rows = np.asarray(reduced.tolist()[:int(rank)], dtype=np.int64)
    pivots = [int(np.flatnonzero(row)[0]) for row in rows]
    return rows, pivots


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spin', type=int, required=True)
    parser.add_argument('--degree', type=int, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=56)
    parser.add_argument('--seconds', type=int, default=2700)
    parser.add_argument('--stagnation-seconds', type=int, default=480)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(parents=True, exist_ok=True)
    assert not (args.output / 'manifest.json').exists()
    start = time.monotonic()
    old_manifest = json.loads((args.source / 'manifest.json').read_text())
    prime = old_manifest['prime']
    weights = tuple(old_manifest['weights'])
    base.setup(weights, (args.spin,), prime)
    data = base.CONTEXT
    folder = args.source / ('s%d' % args.spin)
    old_record = json.loads((folder / 'result.json').read_text())
    assert base.sha(folder / 'payload.pkl') == old_record['payload_sha256']
    payload = pickle.loads((folder / 'payload.pkl').read_bytes())
    selected = [index for index, degree in enumerate(payload['degrees']) if degree == args.degree]
    recipes = [payload['recipes'][index] for index in selected]
    values = list(payload['source_values'][args.degree])
    images = [payload['images'][index] for index in selected]
    assert len(values) == len(recipes) == len(images)
    control_coordinates = 0
    for index in sorted(set([0, len(recipes) - 1])) if recipes else []:
        expression = base.hw_extremal.project({recipes[index]: 1}, (args.spin,) + weights)
        expression = {word: base.Fraction(value).numerator * pow(base.Fraction(value).denominator, -1, prime) % prime
                      for word, value in expression.items() if value}
        assert expression == base.projection.project({recipes[index]: 1}, (args.spin,) + weights, prime)
        q0 = base.native.sp_census.base.action_expr(expression)
        offset = 0
        for source, expected in ((True, values[index]), (False, images[index])):
            frames = {args.degree: data['frames'][args.spin][args.degree]} if source else data['frames'][args.spin - 1]
            offset = 0
            for degree, old_frame in frames.items():
                points = min(2, old_frame.P)
                frame = base.native.Frame(degree, points, 1, prime, k=old_frame.k)
                frame.L = {letter: (odd, matrix[:points].copy()) for letter, (odd, matrix) in old_frame.L.items()}
                parts = {word: value for word, value in (expression if source else q0).items()
                         if sum(3 <= letter <= 5 for trace in word[1] for letter in trace) == degree}
                actual = base.original_evaluate([parts], frame)[0]
                assert np.array_equal(actual, expected[offset:offset + frame.ncols]), ('resume_probe', index, source, degree)
                control_coordinates += frame.ncols
                offset += old_frame.ncols
    dimension = data['dimensions'][args.spin][args.degree]
    source_columns = data['frames'][args.spin][args.degree].ncols
    image_columns = sum(frame.ncols for frame in data['frames'][args.spin - 1].values())
    reduced, pivots = reduce(np.asarray(values, dtype=np.int64).reshape(-1, source_columns), prime)
    assert len(reduced) == len(recipes)
    out, out_pivots = reduce(np.asarray(images, dtype=np.int64).reshape(-1, image_columns), prime)
    initial = len(recipes)
    manifest = dict(group='G2', weights=weights, spin=args.spin, degree=args.degree, prime=prime,
                    exact_degree_dimension=dimension, retained_initial_rows=initial,
                    old_manifest_sha256=base.sha(args.source / 'manifest.json'),
                    old_payload_sha256=base.sha(folder / 'payload.pkl'),
                    original_frames=old_manifest['frames'], source_sha256=base.sha(__file__),
                    base_driver_sha256=base.sha(base.__file__), workers=args.workers,
                    independent_resume_control_coordinates=control_coordinates,
                    scope='One fermion-degree subspace. A combined source-space elimination is required to infer the full Q0 rank or cohomology.')
    base.save(args.output / 'manifest.json', manifest)
    pending, serial, batches = {}, 0, 0
    last, last_gain = time.monotonic(), time.monotonic()

    def write(status):
        path = args.output / 'payload.pkl'
        saved = dict(recipes=recipes, source_values=values, images=images,
                     source_pivots=pivots, outgoing_pivots=out_pivots)
        temporary = path.with_suffix('.tmp')
        temporary.write_bytes(pickle.dumps(saved, protocol=5))
        temporary.replace(path)
        result = dict(status=status, spin=args.spin, degree=args.degree, prime=prime,
                      source_rank=len(recipes), exact_degree_dimension=dimension,
                      degree_basis_complete=len(recipes) == dimension,
                      outgoing_rank_lower_bound=len(out), initial_rows=initial,
                      batches=batches, seconds=time.monotonic() - start,
                      seconds_since_gain=time.monotonic() - last_gain,
                      payload_sha256=base.sha(path), manifest_sha256=base.sha(args.output / 'manifest.json'),
                      job_id=os.environ['SLURM_JOB_ID'], node=os.environ.get('SLURMD_NODENAME'))
        base.save(args.output / 'result.json', result)
        print(json.dumps(result), flush=True)

    workers = min(args.workers, max(1, (dimension - initial + 11) // 12))
    pool = ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('fork'))
    failed = False
    try:
        while pending or len(recipes) < dimension:
            can_submit = (len(recipes) < dimension and time.monotonic() - start < args.seconds
                          and time.monotonic() - last_gain < args.stagnation_seconds)
            while can_submit and len(pending) < workers:
                task = args.spin, args.degree, 2026111801 + 10000019 * args.spin + 100003 * args.degree + 7919 * serial, 12
                pending[pool.submit(base.batch, task)] = serial
                serial += 1
            if not pending:
                break
            done, _ = wait(pending, timeout=1, return_when=FIRST_COMPLETED)
            for future in done:
                pending.pop(future)
                result = future.result()
                batches += 1
                if len(recipes) == dimension:
                    continue
                reduced, pivots, indices = extend_rows(reduced, pivots, result['source_values'], prime)
                assert len(reduced) <= dimension
                if indices:
                    last_gain = time.monotonic()
                    new_images = result['images'][indices]
                    out, out_pivots, _ = extend_rows(out, out_pivots, new_images, prime)
                    recipes.extend(result['recipes'][index] for index in indices)
                    values.extend(result['source_values'][indices])
                    images.extend(new_images)
            if time.monotonic() - last > 60:
                write('running')
                last = time.monotonic()
    except BaseException:
        failed = True
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
        write('failed_checkpointed' if failed else 'complete_physical_degree_basis' if len(recipes) == dimension
              else 'stagnation_guard' if time.monotonic() - last_gain >= args.stagnation_seconds else 'bounded_partial_degree_basis')


if __name__ == '__main__':
    main()
