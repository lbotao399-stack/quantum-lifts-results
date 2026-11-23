"""Rebuild the F4 sixth-power boundary minor and the independent new rows."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import json
from multiprocessing import get_context
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'global_boundaries'))
from adjoint_edge_search import central, old
from albert_boundaries import evaluate
from global_source_boundaries import minor
from mixed_source_search import previous, solve_left
from global_source_scan import square_determinant


def native_old(point, prime):
    return old((point, prime))


def record(weights, values, ids, base_det, prime, fixed=None):
    count = len(ids) - len(weights)
    assert count == 3002 and values.shape[0] == len(ids)
    residual = (values[count:].astype(np.int64) -
                weights.astype(np.int64) @ values[:count].astype(np.int64)) % prime
    if fixed is None:
        rows, cols, complement = minor(residual, prime)
    else:
        rows = [ids.index(i) - count for i in fixed['new_points']]
        cols = fixed['new_columns']
        assert min(rows) >= 0
        complement = square_determinant(residual[np.ix_(rows, cols)], prime)
    assert complement and count + len(rows) <= 3017
    return dict(prime=prime, previous_rank=count, previous_minor=int(base_det),
                new_points=[ids[count+i] for i in rows], new_columns=list(cols),
                complementary_minor=int(complement), rank=count+len(rows),
                determinant=int(base_det*complement%prime))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', required=True, type=int, choices=(1009,1013))
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--fixed', type=Path)
    parser.add_argument('--weights', type=Path, help='Save the verified evaluation functionals here')
    parser.add_argument('--reference', type=Path, help='Save the newly computed reference columns here')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    p = args.prime
    pts = previous()['rows'] + central()['new_rows'] + list(range(7000,7040))
    assert len(pts) == len(set(pts)) == 3042
    with ProcessPoolExecutor(args.workers, mp_context=get_context('fork')) as pool:
        base = np.stack(list(pool.map(partial(native_old, prime=p), pts, chunksize=2)))
        with np.errstate(all='ignore'):
            weights, det = solve_left(base[:3002], base[3002:], p)
        assert int(det) == {1009:833, 1013:615}[p]
        assert np.array_equal(weights.astype(np.int64) @ base[:3002].astype(np.int64) % p, base[3002:])
        print('verified old rows', int(det), flush=True)
        if args.weights:
            np.savez_compressed(args.weights, weights=weights.astype(np.uint16))
        values = np.stack(list(pool.map(partial(evaluate, prime=p), pts, chunksize=2)))
    if args.reference:
        np.savez_compressed(args.reference, points=pts, values=values)
    fixed = json.loads(args.fixed.read_bytes()) if args.fixed else None
    result = record(weights, values, pts, int(det), p, fixed)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
