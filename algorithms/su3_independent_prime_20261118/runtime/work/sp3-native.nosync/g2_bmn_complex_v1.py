"""Sample rationally projected cochains until exact G2 dimensions saturate.

Keep source recipes and Q0 images. Incomplete runs give rank lower bounds;
positive cohomology upper bounds are not certified classes.
"""
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from fractions import Fraction
from math import ceil, comb
from pathlib import Path
import argparse
import hashlib
import json
import multiprocessing
import os
import pickle
import random
import time

import g2_bmn_native_v1 as native
import hw_extremal
import native_bmn_projection_v1 as projection
import numpy as np
from extend_rows_binary64_exact_v1 import extend_rows
from fortuity.bmn.words import contents_of_sector
from sp_projected_basis import evaluate as original_evaluate
from sp3_randomized_quotient_v1 import partitions

CONTEXT = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def setup(weights, spins, prime, points_override=0):
    global CONTEXT
    native.configure(prime)
    dimensions, frames, metadata, choices = {}, {}, {}, {}
    for spin in range(max(spins) + 1):
        dimensions[spin] = {int(degree): int(dimension) for degree, dimension in
                            native.sp_census.base.exact_hw((spin,) + weights, 0).items() if dimension}
        frames[spin], metadata[spin] = {}, {}
        for degree, dimension in sorted(dimensions[spin].items()):
            k = degree + 2 if degree else 0
            points = points_override or ceil((dimension + 24) / comb(k, degree))
            seed = 2026111787 + 100003 * spin + 7919 * degree
            frames[spin][degree] = native.Frame(degree, points, seed, prime, k=k)
            metadata[spin][degree] = dict(degree=degree, points=points, k=k, seed=seed,
                                          columns=frames[spin][degree].ncols)
        if spin not in spins:
            continue
        choices[spin] = defaultdict(list)
        for content in contents_of_sector((spin,) + weights):
            degree = sum(content[3:6])
            if degree in dimensions[spin]:
                letters = tuple(index for index, count in enumerate(content) for _ in range(count))
                choices[spin][degree].append((letters, tuple(partitions(len(letters)))))
    CONTEXT = dict(weights=weights, spins=spins, prime=prime, dimensions=dimensions,
                   frames=frames, frame_metadata=metadata, choices=choices)


def batch(task):
    spin, degree, seed, count = task
    data = CONTEXT
    prime = data['prime']
    rng = random.Random(seed)
    recipes, expressions, seen = [], [], set()
    attempts = 0
    start = time.monotonic()
    while len(recipes) < count and attempts < count * 120:
        candidates = data['choices'][spin][degree]
        letters, shapes = candidates[(attempts + seed) % len(candidates)]
        attempts += 1
        letters = list(letters)
        rng.shuffle(letters)
        shape = rng.choice(shapes)
        traces, offset = [], 0
        for length in shape:
            traces.append(tuple(letters[offset:offset + length]))
            offset += length
        sign, word = native.sp_census.base.canon((), tuple(traces))
        if not sign or word in seen:
            continue
        seen.add(word)
        expression = projection.project({word: 1}, (spin,) + data['weights'], prime)
        if expression:
            recipes.append(word)
            expressions.append(expression)
    projected = time.monotonic()
    values = native.evaluate(expressions, data['frames'][spin][degree])
    images = np.concatenate([native.evaluate(expressions, frame, derivative=True)
                             for frame in data['frames'][spin - 1].values()], axis=1)
    return dict(spin=spin, degree=degree, recipes=recipes, source_values=values, images=images,
                attempts=attempts, projected_terms=sum(map(len, expressions)),
                projection_seconds=projected - start, evaluation_seconds=time.monotonic() - projected,
                seconds=time.monotonic() - start)


def control(output):
    records = []
    prime = CONTEXT['prime']
    for spin in CONTEXT['spins']:
        for degree in sorted(CONTEXT['dimensions'][spin]):
            result = batch((spin, degree, 91339 + 1009 * spin + 71 * degree, 2))
            source, target = [], []
            for word in result['recipes']:
                expression = hw_extremal.project({word: 1}, (spin,) + CONTEXT['weights'])
                reduced = {word: Fraction(value).numerator * pow(Fraction(value).denominator, -1, prime) % prime
                           for word, value in expression.items() if value}
                assert reduced == projection.project({word: 1}, (spin,) + CONTEXT['weights'], prime)
                assert not hw_extremal.action(expression, 0, 1)
                assert not hw_extremal.action(expression, 1, 2)
                image = native.sp_census.base.action_expr(expression)
                assert not native.sp_census.base.action_expr(image)
                image = {word: Fraction(value).numerator * pow(Fraction(value).denominator, -1, prime) % prime
                         for word, value in image.items() if value}
                source.append(original_evaluate([reduced], CONTEXT['frames'][spin][degree])[0])
                blocks = []
                for target_degree, frame in CONTEXT['frames'][spin - 1].items():
                    part = {word: value for word, value in image.items()
                            if sum(3 <= letter <= 5 for trace in word[1] for letter in trace) == target_degree}
                    blocks.append(original_evaluate([part], frame)[0])
                target.append(np.concatenate(blocks))
            if source:
                assert np.array_equal(np.asarray(source), result['source_values']), (spin, degree, 'source')
                assert np.array_equal(np.asarray(target), result['images']), (spin, degree, 'Q0')
            records.append(dict(spin=spin, degree=degree, rows=len(source),
                                source_coordinates=result['source_values'].size, Q0_coordinates=result['images'].size,
                                nonzero_Q0_coordinates=int(np.count_nonzero(result['images']))))
    assert sum(row['nonzero_Q0_coordinates'] for row in records) > 0
    result = dict(status='passed_rational_projection_and_independent_symbolic_Q0_controls',
                  weights=CONTEXT['weights'], prime=prime, records=records,
                  lie_algebra=native.g2_basis()[1], driver_sha256=sha(__file__), backend_sha256=sha(native.__file__))
    save(output / 'control.json', result)
    print(json.dumps(result), flush=True)


def run(args):
    assert os.environ.get('SLURM_JOB_ID'), 'Production cochains run on the authorized cluster'
    start = time.monotonic()
    data, prime = CONTEXT, args.prime
    states, pending, serial = {}, {}, 0
    for spin in data['spins']:
        dimensions = data['dimensions'][spin]
        columns = sum(frame.ncols for frame in data['frames'][spin - 1].values())
        states[spin] = dict(source={degree: (np.zeros((0, data['frames'][spin][degree].ncols), dtype=np.int64), [])
                                   for degree in dimensions},
                            outgoing=(np.zeros((0, columns), dtype=np.int64), []),
                            recipes=[], degrees=[], images=[], source_values={degree: [] for degree in dimensions},
                            selected_outgoing=[], batches=0, attempts=0, projected_terms=0,
                            projection_seconds=0.0, evaluation_seconds=0.0)
    manifest = dict(group='G2', weights=data['weights'], prime=prime, dimensions=data['dimensions'],
                    frames=data['frame_metadata'], workers=args.workers, wall_seconds=args.seconds,
                    sources={Path(path).name: sha(path) for path in
                             (__file__, native.__file__, projection.__file__,
                              native.HERE / 'bmn_derivation_rank_v1.cpp', native.HERE / 'trace_frame_bmn_rank_v1.cpp',
                              native.HERE / ('bmn_derivation_rank_v1_d7_p%d.dylib' % prime))},
                    rational_bridge='Each retained source is the pinned rational extremal projection of its saved integer trace word. Its reduction agrees modulo the prime. Independent evaluation minors certify rational ranks; character dimensions certify completeness only upon saturation.')
    save(args.output / 'manifest.json', manifest)

    def write(status):
        ranks = {spin: len(state['outgoing'][0]) for spin, state in states.items()}
        records = []
        for spin, state in states.items():
            folder = args.output / ('s%d' % spin)
            folder.mkdir(exist_ok=True)
            complete = all(len(state['source'][degree][0]) == dimension for degree, dimension in data['dimensions'][spin].items())
            payload = {key: state[key] for key in ('recipes', 'degrees', 'images', 'source_values', 'selected_outgoing')}
            payload['source_pivots'] = {degree: values[1] for degree, values in state['source'].items()}
            payload['outgoing_pivots'] = state['outgoing'][1]
            path = folder / 'payload.pkl'
            temporary = path.with_suffix('.tmp')
            temporary.write_bytes(pickle.dumps(payload, protocol=5))
            temporary.replace(path)
            dimension = sum(data['dimensions'][spin].values())
            incoming = ranks.get(spin + 1, 0)
            upper = dimension - ranks[spin] - incoming
            assert upper >= 0
            record = dict(sector=(spin,) + data['weights'], dimension=dimension,
                          retained_source_rows=len(state['recipes']), complete_physical_basis=complete,
                          source_ranks={degree: len(value[0]) for degree, value in state['source'].items()},
                          outgoing_rank_lower_bound=ranks[spin], incoming_rank_lower_bound=incoming,
                          cohomology_upper_bound=upper, rational_zero_cohomology_certified=upper == 0,
                          payload_sha256=sha(path), **{key: state[key] for key in
                          ('batches', 'attempts', 'projected_terms', 'projection_seconds', 'evaluation_seconds')})
            save(folder / 'result.json', record)
            records.append(record)
        summary = dict(status=status, group='G2', weights=data['weights'], prime=prime, records=records,
                       scalar_cohomology_upper_bound=sum(data['dimensions'][0].values()) - ranks.get(1, 0),
                       seconds=time.monotonic() - start, job_id=os.environ['SLURM_JOB_ID'],
                       scope='Exact character dimensions and rational rank lower bounds. Complete physical bases are recorded only at dimension saturation. Positive remainders are not certified classes or quantum obstructions.')
        save(args.output / 'result.json', summary)
        print(json.dumps(dict(status=status, retained=[row['retained_source_rows'] for row in records],
                              ranks=ranks, seconds=summary['seconds'])), flush=True)

    def next_task():
        counts = Counter((task[0], task[1]) for task in pending.values())
        choices = []
        for spin, state in states.items():
            for degree, dimension in data['dimensions'][spin].items():
                missing = dimension - len(state['source'][degree][0])
                if missing:
                    # Balance sectors first, then their unresolved degrees.
                    weight = missing / dimension / (1 + sum(value for (s, _), value in counts.items() if s == spin))
                    weight /= 1 + counts[spin, degree]
                    choices.append((weight, -spin, -degree, spin, degree))
        return max(choices)[-2:] if choices else None

    pool = ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('fork'))
    last, failed = time.monotonic(), False
    try:
        while pending or time.monotonic() - start < args.seconds:
            while len(pending) < args.workers and time.monotonic() - start < args.seconds:
                choice = next_task()
                if choice is None:
                    break
                spin, degree = choice
                task = spin, degree, 2026111793 + 7919 * serial, args.batch_size
                pending[pool.submit(batch, task)] = task
                serial += 1
            if not pending:
                break
            done, _ = wait(pending, timeout=1, return_when=FIRST_COMPLETED)
            for future in done:
                pending.pop(future)
                result = future.result()
                spin, degree = result['spin'], result['degree']
                state = states[spin]
                state['batches'] += 1
                for key in ('attempts', 'projected_terms', 'projection_seconds', 'evaluation_seconds'):
                    state[key] += result[key]
                reduced, pivots = state['source'][degree]
                reduced, pivots, indices = extend_rows(reduced, pivots, result['source_values'], prime)
                assert len(reduced) <= data['dimensions'][spin][degree]
                state['source'][degree] = reduced, pivots
                if indices:
                    old = len(state['recipes'])
                    images = result['images'][indices]
                    out, out_pivots = state['outgoing']
                    out, out_pivots, independent = extend_rows(out, out_pivots, images, prime)
                    state['outgoing'] = out, out_pivots
                    state['selected_outgoing'].extend(old + index for index in independent)
                    state['recipes'].extend(result['recipes'][index] for index in indices)
                    state['degrees'].extend([degree] * len(indices))
                    state['images'].extend(images)
                    state['source_values'][degree].extend(result['source_values'][indices])
            if time.monotonic() - last > 45:
                write('running')
                last = time.monotonic()
    except BaseException:
        failed = True
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
        complete = all(len(state['source'][degree][0]) == dimension for spin, state in states.items()
                       for degree, dimension in data['dimensions'][spin].items())
        write('failed_checkpointed' if failed else 'complete_physical_cochain_bases' if complete else 'bounded_partial_physical_bases')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--weights', default='5,5,5')
    parser.add_argument('--spins', default='1,2,3,4,5,6,7,8,9')
    parser.add_argument('--prime', type=int, default=1000003)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=36)
    parser.add_argument('--batch-size', type=int, default=12)
    parser.add_argument('--seconds', type=int, default=900)
    parser.add_argument('--control', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    assert not (args.output / 'manifest.json').exists()
    setup(tuple(map(int, args.weights.split(','))), tuple(map(int, args.spins.split(','))), args.prime,
          points_override=2 if args.control else 0)
    if args.control:
        control(args.output)
    else:
        run(args)


if __name__ == '__main__':
    main()
