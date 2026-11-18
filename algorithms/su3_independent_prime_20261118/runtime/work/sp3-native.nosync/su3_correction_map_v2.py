"""Checkpointed blocks of the coupled first-correction equation.

These are evaluation matrices, not yet injectivity or rational-solution
certificates. All structural Q0 blocks, including zero-RHS blocks, are kept.
"""
import sp_environment
import su_trace_adapter_v1 as algebra
from native_su_projection_v1 import project
from native_su_parallel_v1 import evaluate_profile
from g2_full_classical_screen_v1 import sector
from fortuity.full.evalg import Frame
from fortuity.full.letters import letters_up_to
from pathlib import Path
from math import ceil, comb
import argparse
import hashlib
import json
import pickle
import time
import numpy as np

HERE = Path(__file__).resolve().parent
PRIME = 1000033


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path = Path(path)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data, indent=2)+'\n')
    temp.replace(path)


def target_frame(target, seed=2026111791):
    degree = target['fermions']
    k = degree + 2
    points = ceil((target['dimension']+32)/comb(k, degree))
    alphabet = sorted(L for L in letters_up_to(2, 2) if sum(L[:2]) <= 2)
    actual_seed = seed + 1009*target['index']
    frame = Frame(3, points, k, [degree], alphabet, seed=actual_seed, p=PRIME)
    return frame, dict(N=3, points=points, k=k, degrees=[degree], alphabet=alphabet,
                      seed=actual_seed, prime=PRIME, columns=frame.ncols)


def run(index, basis, output, start_row=0, stop_row=0, batch=8, seconds=600, workers=8):
    started = time.monotonic()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plan = json.loads((HERE/'su3_correction_profiles_v1.json').read_text())
    source = plan['sources'][index]
    source_path = Path(basis)/('profile%02d'%index)/'independent_seeds.pkl'
    seeds = pickle.loads(source_path.read_bytes())
    available = len(seeds)
    stop_row = stop_row or available
    if stop_row > available:
        result = dict(status='deferred_missing_source_rows', source=index, start_row=start_row,
                      stop_row=stop_row, available_rows=available, correction_solved=False)
        save(output/'deferred.json', result)
        return result
    assert 0 <= start_row < stop_row <= available
    seeds = seeds[start_row:stop_row]
    recipes = pickle.dumps(seeds, protocol=4)
    recipe_path = output/'source_recipes.pkl'
    if recipe_path.exists():
        assert recipe_path.read_bytes() == recipes
    else:
        recipe_path.write_bytes(recipes)
    frames = {j: target_frame(plan['targets'][j]) for j in source['targets']}
    manifest = dict(source=source, prime=PRIME, mu=plan['mu'], length=9, parity=-1,
                    source_recipes_sha256=sha(recipe_path), start_row=start_row, stop_row=stop_row, rows=len(seeds), complete_source_dimension=source['dimension'],
                    targets={j: description for j, (_, description) in frames.items()}, batch=batch,
                    driver_sha256=sha(__file__), parallel_driver_sha256=sha(HERE/'native_su_parallel_v1.py'),
                    plan_sha256=sha(HERE/'su3_correction_profiles_v1.json'))
    manifest = json.loads(json.dumps(manifest))
    manifest_path = output/'manifest.json'
    if manifest_path.exists():
        assert json.loads(manifest_path.read_text()) == manifest
    else:
        save(manifest_path, manifest)
    records = []
    timed_out = False
    for start in range(start_row, stop_row, batch):
        stop = min(start+batch, stop_row)
        missing = [j for j in source['targets'] if not (output/('t%02d_%05d_%05d.json'%(j,start,stop))).exists()]
        if missing:
            before = time.monotonic()
            expressions = [project(algebra.parity_project({word:1}, -1), sector(tuple(plan['mu']), 9), PRIME)
                           for word in seeds[start-start_row:stop-start_row]]
            projection_seconds = time.monotonic()-before
        for j in source['targets']:
            name = 't%02d_%05d_%05d'%(j,start,stop)
            marker, matrix_path = output/(name+'.json'), output/(name+'.npy')
            if marker.exists():
                record = json.loads(marker.read_text())
                assert sha(matrix_path) == record['sha256']
                records.append(record)
                continue
            frame = frames[j][0]
            before = time.monotonic()
            matrix = evaluate_profile(expressions, frame, source['profile'], plan['targets'][j]['profile'], workers)
            elapsed = time.monotonic()-before
            assert matrix.shape == (stop-start, frame.ncols)
            temporary = matrix_path.with_suffix('.tmp.npy')
            np.save(temporary, matrix.astype(np.uint32))
            temporary.replace(matrix_path)
            record = dict(source=index, target=j, start=start, stop=stop, columns=frame.ncols,
                          projection_seconds=projection_seconds, evaluation_seconds=elapsed,
                          projected_terms=sum(map(len, expressions)), nonzero_entries=int(np.count_nonzero(matrix)),
                          workers=workers, sha256=sha(matrix_path))
            save(marker, record)
            records.append(record)
            print(json.dumps(record), flush=True)
            save(output/'progress.json', dict(source=index, completed_blocks=len(records),
                 required_blocks=ceil(len(seeds)/batch)*len(source['targets']), latest=record,
                 seconds=time.monotonic()-started))
        if time.monotonic()-started >= seconds:
            timed_out = True
            break
    required = ceil(len(seeds)/batch)*len(source['targets'])
    result = dict(status='complete_selected_evaluation_blocks' if len(records)==required else 'partial_evaluation_blocks',
                  source=index, start_row=start_row, stop_row=stop_row, input_rows=len(seeds), complete_source_dimension=source['dimension'],
                  completed_blocks=len(records), required_blocks=required, records=records,
                  seconds=time.monotonic()-started, stopped_by='time_cap' if timed_out and len(records)<required else 'selected_rows_done',
                  all_structural_targets=source['targets'], correction_solved=False,
                  scope='Evaluation blocks only; full target-frame injectivity and the coupled solve remain to be certified.')
    save(output/'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--basis', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--start-row', type=int, default=0)
    parser.add_argument('--stop-row', type=int, default=0)
    parser.add_argument('--batch', type=int, default=8)
    parser.add_argument('--seconds', type=int, default=600)
    parser.add_argument('--workers', type=int, default=8)
    print(json.dumps(run(**vars(parser.parse_args()))), flush=True)
