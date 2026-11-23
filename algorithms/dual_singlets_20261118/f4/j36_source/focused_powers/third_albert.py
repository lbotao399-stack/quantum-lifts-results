"""Classical maps on the five alternating triples of the sixth glueball's third source."""

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import combinations, permutations
from multiprocessing import get_context
from pathlib import Path
import argparse
import json

import numpy as np

from albert_boundaries import SIGNS, action, apply, balanced, colour_tensors, features, native26
from sandwich import HERE, WIDTH, first, minor, run_new, run_old, solve_left


def triple(tensor, scalars, marked, axes, prime):
    total = np.zeros_like(tensor)
    if marked is None:
        for i, j, k in permutations(range(3)):
            start = balanced(apply(scalars[i], tensor, (axes[0],)), prime)
            total += SIGNS[i, j, k] * apply(scalars[k], apply(scalars[j], start, (axes[1],)), (axes[2],))
    else:
        for pos, (i, j) in enumerate(combinations(range(3), 2)):
            k = 3 - i - j
            start = apply(marked[pos], tensor, (axes[0],))
            total += SIGNS[i, j, k] * apply(scalars[k], start, (axes[1],))
    return balanced(total, prime)


def arms(data, side, order, prime):
    low, high, _, _ = colour_tensors()
    corner = (low, high)[side]
    x, y, comm, jac = data
    if side == 0:
        x, y, comm = [v.transpose(0, 2, 1) for v in (x, y, comm)]
        jac = jac.T
    xs, ys, ks = [[action(a)[0] for a in v] for v in (x, y, comm)]
    mark_axes = [(0, 2), (1, 3)][side]
    scalar_axes = [((0, 1, 3),), ((2, 3, 1), (0, 2, 1))][side]
    starts = {}
    if order == 1 and side == 1:
        starts["plain"] = triple(corner, xs, None, (0, 2, 3), prime)
    else:
        starts["split"] = triple(corner, xs, ks, mark_axes, prime)
    if order == 2:
        starts["mark"] = triple(corner, xs, ys, mark_axes, prime)
        if side == 0:
            starts["jacobi"] = balanced(apply(action(jac)[0], corner, (1,)), prime)
        else:
            starts["plain"] = triple(corner, xs, None, (0, 2, 3), prime)
    for key in starts:
        for axes in scalar_axes:
            starts[key] = triple(starts[key], xs, None, axes, prime)
    return starts


def evaluate(number, prime, order, closed=False):
    x = (np.random.default_rng(741411 + number).integers(-3, 4, size=(6, 52))
         @ native26().reshape(52, -1)).reshape(6, 26, 26)
    def bracket(a, b):
        return balanced(a @ b - b @ a, prime)
    comm = np.stack([bracket(x[i], x[j]) for i, j in combinations(range(3), 2)])
    y = comm if closed else x[3:]
    jac = balanced(bracket(x[0], y[2]) - bracket(x[1], y[1]) + bracket(x[2], y[0]), prime)
    _, _, exchange, transposed = colour_tensors()
    all_arms = []
    for side, transfer in enumerate((transposed, exchange)):
        initial = arms((x[:3], y, comm, jac), side, order, prime)
        all_arms.append(dict(zip(initial, features(initial.values(), transfer, prime))))
    left, right = all_arms
    with np.errstate(all="ignore"):
        if order == 1:
            scores = [left["split"] @ right["plain"].T]
        else:
            scores = [left["split"] @ right["mark"].T - left["mark"] @ right["split"].T,
                      left["jacobi"] @ right["plain"].T]
    scores = np.stack(scores)
    assert np.all(np.isfinite(scores)) and np.array_equal(scores, scores.astype(np.int64))
    return np.asarray(scores.astype(np.int64) % prime, dtype=np.uint16).ravel()


def old_row(number, prime, order, selected):
    base = run_old((3, order, number, prime))
    batches = {c: run_new((3, order, number, prime, c)) for c in set(j // WIDTH for j in selected)}
    return np.r_[base, [batches[j // WIDTH][j % WIDTH] for j in selected]]


def collect(order, prime, workers, certificate):
    prev = json.loads((HERE / f"sandwich-3-{order}.json").read_bytes())
    rows = first(3, order)[0] + prev["new_rows"]
    n = len(rows)
    assert n == prev["rank"]
    test_rows = certificate["new_rows"] if certificate else list(range(7900, 7920))
    ids = rows + test_rows
    with ProcessPoolExecutor(workers, mp_context=get_context("fork")) as pool:
        base = np.array(list(pool.map(partial(old_row, prime=prime, order=order,
                                               selected=prev["new_columns"]), ids, chunksize=2)))
        with np.errstate(all="ignore"):
            weights, det = solve_left(base[:n], base[n:], prime)
        expected = prev["determinant"] if prime == 1009 else json.loads(
            (HERE / f"sandwich-replay-3-{order}.json").read_bytes())["determinant"]
        assert det == expected
        print("verified base", n, det, flush=True)
        arr = np.array(list(pool.map(partial(evaluate, prime=prime, order=order), ids, chunksize=2)))
    residual = (arr[n:].astype(np.int64) - weights.astype(np.int64) @ arr[:n].astype(np.int64)) % prime
    if certificate:
        from global_source_scan import square_determinant
        cols = certificate["new_columns"]
        d = square_determinant(residual[:, cols], prime) if cols else 1
        assert d
        r = list(range(len(test_rows)))
    else:
        r, cols, d = minor(residual, prime)
    return dict(prime=prime, order=order, base_rank=n, rank=n+len(r), base_minor=int(det),
                new_rows=[test_rows[i] for i in r], new_columns=cols, complement_minor=int(d),
                determinant=int(d*det%prime))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("controls", "search", "replay"))
    parser.add_argument("--order", type=int, choices=(1, 2), required=True)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--certificate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.action == "controls":
        for i in (12, 219):
            q = evaluate(i, args.prime, args.order)
            if args.order == 2:
                assert not evaluate(i, args.prime, args.order, True).any()
            assert q.any()
            print(i, len(q), int(np.count_nonzero(q)))
    else:
        record = json.loads(args.certificate.read_bytes()) if args.action == "replay" else None
        result = collect(args.order, args.prime, args.workers, record)
        if args.output:
            args.output.write_text(json.dumps(result, indent=2)+"\n")
        print(result, flush=True)
