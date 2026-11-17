"""Merge disjoint source-degree continuations and test generic G2 products.

All matrix rows use the original frames. Degree slices replace, rather than
add to, the original source slices. Quantum differentials are not computed.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import multiprocessing
import os
import pickle
import time
import zipfile

import g2_bmn_complex_v1 as base
import numpy as np
from flint import nmod_mat
from sp_graviton_singlets import product_values

FAMILIES = {}


def select_minor(values, prime):
    if not len(values) or not values.shape[1]:
        return [], [], 1
    reduced, rank = nmod_mat(values.T.tolist(), prime).rref()
    rows = [next(index for index, value in enumerate(row) if value)
            for row in reduced.tolist()[:int(rank)]]
    if not rows:
        return [], [], 1
    reduced, rank2 = nmod_mat(values[rows].tolist(), prime).rref()
    assert rank2 == rank
    columns = [next(index for index, value in enumerate(row) if value)
               for row in reduced.tolist()[:int(rank)]]
    determinant = int(nmod_mat(values[np.ix_(rows, columns)].tolist(), prime).det())
    assert determinant
    return rows, columns, determinant


def product_block(task):
    spin, degree, first, last = task
    family = FAMILIES[spin]
    original = base.CONTEXT['frames'][spin][degree]
    frame = base.native.Frame(degree, last - first, 1, original.p, k=original.k)
    frame.L = {letter: (odd, matrix[first:last].copy())
               for letter, (odd, matrix) in original.L.items()}
    kernel = family['kernel']
    needed = sorted(set(map(int, kernel.indices)))
    assert int(max(np.diff(kernel.indptr), default=0)) * (original.p - 1) ** 2 < 2 ** 63
    raw = product_values(family['generators'], [family['products'][i] for i in needed], frame, degree)
    values = np.asarray(kernel[:, needed] @ raw) % original.p
    return task, values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--degree-runs', type=Path, required=True)
    parser.add_argument('--family-run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=48)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    assert not args.output.exists()
    args.output.mkdir(parents=True)
    start = time.monotonic()
    manifest = json.loads((args.source / 'manifest.json').read_text())
    prime = manifest['prime']
    base.setup(tuple(manifest['weights']), tuple(range(1, 10)), prime)
    assert json.loads(json.dumps(base.CONTEXT['dimensions'])) == manifest['dimensions']
    assert json.loads(json.dumps(base.CONTEXT['frame_metadata'])) == manifest['frames']
    cochains, ranks, selections, records = {}, {0: 0, 10: 0}, {}, []
    evidence = args.output / 'evidence'
    evidence.mkdir()
    base.save(evidence / 'frame_manifest.json', manifest)
    for spin in range(1, 10):
        folder = args.source / ('s%d' % spin)
        record = json.loads((folder / 'result.json').read_text())
        assert base.sha(folder / 'payload.pkl') == record['payload_sha256']
        old = pickle.loads((folder / 'payload.pkl').read_bytes())
        recipes, degrees, images, replacements = [], [], [], []
        for degree in sorted(base.CONTEXT['dimensions'][spin]):
            indices = [i for i, value in enumerate(old['degrees']) if value == degree]
            new = args.degree_runs / ('s%d_d%d' % (spin, degree))
            if new.exists():
                state = json.loads((new / 'result.json').read_text())
                assert state['status'] != 'running'
                assert base.sha(new / 'payload.pkl') == state['payload_sha256']
                extension = pickle.loads((new / 'payload.pkl').read_bytes())
                assert extension['recipes'][:len(indices)] == [old['recipes'][i] for i in indices]
                assert np.array_equal(np.asarray(extension['images'][:len(indices)]),
                                      np.asarray([old['images'][i] for i in indices]))
                recipes.extend(extension['recipes'])
                degrees.extend([degree] * len(extension['recipes']))
                images.extend(extension['images'])
                replacements.append(dict(degree=degree, result=state, manifest_sha256=base.sha(new / 'manifest.json')))
            else:
                recipes.extend(old['recipes'][i] for i in indices)
                degrees.extend([degree] * len(indices))
                images.extend(old['images'][i] for i in indices)
        images = np.asarray(images, dtype=np.int64)
        rows, columns, determinant = select_minor(images, prime)
        ranks[spin] = len(rows)
        assert ranks[spin] >= record['outgoing_rank_lower_bound']
        cochains[spin] = dict(recipes=recipes, degrees=degrees, images=images)
        selections[spin] = rows
        # Two-point production probes permit independent symbolic Q0 replay.
        probes, offset = [], 0
        selected_probes = sorted(set(rows[:1] + rows[len(rows)//2:len(rows)//2+1] + rows[-1:]))
        for degree, frame in base.CONTEXT['frames'][spin - 1].items():
            points = min(2, frame.P)
            width = points * frame.ncols // frame.P
            probes.append(dict(degree=degree, k=frame.k, points=points,
                letters={letter: (odd, matrix[:points].copy()) for letter, (odd, matrix) in frame.L.items()},
                values=images[np.ix_(selected_probes, range(offset, offset + width))].copy()))
            offset += frame.ncols
        certificate = dict(spin=spin, prime=prime, recipes=[recipes[i] for i in rows],
            minor=images[np.ix_(rows, columns)].copy(), columns=columns, determinant=determinant,
            probe_recipes=[recipes[i] for i in selected_probes], probes=probes)
        (evidence / ('s%d.pkl' % spin)).write_bytes(pickle.dumps(certificate, protocol=5))
        records.append(dict(spin=spin, dimension=sum(base.CONTEXT['dimensions'][spin].values()),
            retained_source_rows=len(recipes), outgoing_rank_lower_bound=ranks[spin], replacements=replacements))
        print(json.dumps(dict(stage='merged_rank', spin=spin, rank=ranks[spin], rows=len(recipes))), flush=True)
    dimensions = {spin: sum(value.values()) for spin, value in base.CONTEXT['dimensions'].items()}
    bounds = []
    for spin in range(10):
        upper = dimensions[spin] - ranks[spin] - ranks[spin + 1]
        assert upper >= 0
        bounds.append(dict(spin=spin, dimension=dimensions[spin],
            outgoing_rank_lower_bound=ranks[spin], incoming_rank_lower_bound=ranks[spin + 1],
            cohomology_upper_bound=upper, rational_zero_cohomology_certified=upper == 0))
    base.save(args.output / 'merged_ranks.json', dict(records=records, bounds=bounds))
    blocks, tasks = {}, []
    closure = set()
    for spin in (1, 3, 5, 7):
        path = args.family_run / ('s%d_5_5_5' % spin) / 'family.pkl'
        family = pickle.loads(path.read_bytes())
        FAMILIES[spin] = family
        for generator in family['generators']:
            fingerprint = repr(sorted(generator['expr'].items()))
            if fingerprint not in closure:
                assert not base.native.sp_census.base.action_expr(generator['expr'])
                closure.add(fingerprint)
        blocks[spin] = {}
        for degree, frame in base.CONTEXT['frames'][spin].items():
            blocks[spin][degree] = np.zeros((family['kernel'].shape[0], frame.ncols), dtype=np.int64)
            if family['kernel'].shape[0]:
                step = max(1, (frame.P + 7) // 8)
                tasks.extend((spin, degree, first, min(frame.P, first + step))
                             for first in range(0, frame.P, step))
        (evidence / ('family_s%d.pkl' % spin)).write_bytes(path.read_bytes())
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=multiprocessing.get_context('fork')) as pool:
        for future in as_completed([pool.submit(product_block, task) for task in tasks]):
            (spin, degree, first, last), values = future.result()
            frame = base.CONTEXT['frames'][spin][degree]
            width = frame.ncols // frame.P
            blocks[spin][degree][:, first * width:last * width] = values
    product_records = []
    for spin in (1, 3, 5, 7):
        values = np.concatenate(list(blocks[spin].values()), axis=1)
        incoming = cochains[spin + 1]['images'][selections[spin + 1]]
        combined = np.concatenate([incoming, values], axis=0)
        rows, columns, determinant = select_minor(combined, prime)
        incoming_upper = dimensions[spin + 1] - ranks.get(spin + 2, 0)
        lower = max(0, len(rows) - incoming_upper)
        assert lower <= bounds[spin]['cohomology_upper_bound']
        row = dict(spin=spin, formal_highest_weight_rows=len(values),
            incoming_rank_lower_bound=ranks[spin + 1], incoming_rank_upper_bound=incoming_upper,
            combined_rank_lower_bound=len(rows), known_product_cohomology_lower_bound=lower,
            cohomology_upper_bound=bounds[spin]['cohomology_upper_bound'],
            known_products_exhaust_cohomology=lower == bounds[spin]['cohomology_upper_bound'])
        product_records.append(row)
        (evidence / ('products_s%d.pkl' % spin)).write_bytes(pickle.dumps(dict(
            values=values, minor=combined[np.ix_(rows, columns)], rows=rows, columns=columns,
            determinant=determinant, incoming_recipes=[cochains[spin+1]['recipes'][i] for i in selections[spin+1]],
            result=row), protocol=5))
        print(json.dumps(dict(stage='generic_products', **row)), flush=True)
    summary = dict(status='completed_merged_degree_ranks_and_generic_products', group='G2', weights=[5,5,5],
        prime=prime, bounds=bounds, products=product_records, exact_closed_generators=len(closure),
        seconds=time.monotonic()-start, driver_sha256=base.sha(__file__),
        scope='Rational rank lower bounds and exact dimension upper bounds. Known-product quotient tests use an independent upper bound on incoming rank. No quantum differential or new-generator claim.')
    base.save(evidence / 'result.json', summary)
    hashes = {str(p.relative_to(evidence)): base.sha(p) for p in evidence.rglob('*') if p.is_file()}
    base.save(evidence / 'SHA256SUMS.json', hashes)
    archive = args.output / 'evidence.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
        for path in sorted(evidence.iterdir()):
            z.write(path, path.name)
    summary['archive'] = dict(path=str(archive), bytes=archive.stat().st_size, sha256=base.sha(archive))
    base.save(args.output / 'result.json', summary)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
