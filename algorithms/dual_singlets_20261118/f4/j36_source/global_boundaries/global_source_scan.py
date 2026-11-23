"""Combine the certified outer minors with direct Jacobi and Koszul predecessors."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from multiprocessing import get_context
import json
from pathlib import Path
import numpy as np
from global_source_boundaries import HERE, colour, columns, minor, write
from global_source_jacobi import JacobiEvaluation

# Each additional test uses a fixed reproducible subset of actual precursors.
FAMILIES = (
    ("in", 0),
    ("in", 1),
    ("in", 2),
    ("out", 2),
    ("in", 3),
    ("out", 3),
    ("in", 4),
    ("in", 5),
    ("in", 6),
    ("out", 6),
    ("in", 8),
    ("out", 8),
    ("in", 7),
)
WIDTH = 850


@cache
def previous():
    saved = json.loads((HERE / "global-2496.json").read_bytes())
    assert saved["rank_lower_bound"] == 2398
    choices = np.asarray(saved["columns"])
    width = saved["width"]
    order = []
    for j, case in enumerate(saved["cases"]):
        positions = np.flatnonzero(choices // width == j)
        order.append(
            (case, positions, columns(case, width)[choices[positions] % width])
        )
    return order


@cache
def sample_columns(kind, case):
    size = 4 * (42 * 65 if kind == "in" else 64 * 48)
    return np.sort(
        np.random.default_rng(414039 + 100 * case + int(kind == "out")).choice(
            size, WIDTH, replace=False
        )
    )


def sample(task):
    number, p = task
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        e = JacobiEvaluation(14337 + number, p)
        old = np.empty(2398, dtype=np.uint16)
        for case, positions, coords in previous():
            old[positions] = e.block(case)[coords]
        result = [old]
        for kind, c in FAMILIES:
            values = e.internal(c) if kind == "in" else e.external(c)
            result.append(values[sample_columns(kind, c)])
        return np.concatenate(result)


def square_determinant(matrix, p):
    # Independent row-pivoting block LU in fixed column order.
    a = matrix.astype(float)
    n = len(a)
    det = 1
    for start in range(0, n, 32):
        end = min(n, start + 32)
        width = end - start
        left = a[start:, start:end].astype(np.int64)
        top = []
        for k in range(width):
            choices = np.flatnonzero(left[k:, k])
            if not len(choices):
                return 0
            i = k + int(choices[0])
            if i != k:
                left[[k, i]] = left[[i, k]]
                a[[start + k, start + i]] = a[[start + i, start + k]]
                det = -det
                for row in top:
                    row[k], row[i] = row[i], row[k]
            lead = int(left[k, k])
            det = det * lead % p
            factors = left[:, k] * pow(lead, -1, p) % p
            factors[: k + 1] = 0
            left = (left - factors[:, None] * left[k]) % p
            top.append(factors)
        u = a[start:end, end:].copy()
        for k, row in enumerate(top):
            if k + 1 < width:
                u[k + 1 :] = (u[k + 1 :] - row[k + 1 : width, None] * u[k]) % p
        if end < n:
            multipliers = np.stack([v[width:] for v in top], axis=1).astype(float)
            a[end:, end:] = (a[end:, end:] - multipliers @ u) % p
    return det


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("search", "replay"))
    parser.add_argument("--certificate", type=Path)
    parser.add_argument("--out", type=Path, default=Path("."))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--points", type=int, default=3088)
    args = parser.parse_args()
    assert __debug__ and args.workers >= 1
    prime = args.prime
    colour(prime)
    saved = None
    if args.action == "replay":
        saved = json.loads(args.certificate.read_bytes())
        numbers = saved["rows"]
        cols = saved["columns"]
        width = len(cols)
    else:
        numbers = range(args.points)
        width = 2398 + WIDTH * len(FAMILIES)
    a = np.empty((len(numbers), width), dtype=np.uint16)
    limits = {2528, 2816, len(numbers)}
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        for i, row in enumerate(
            pool.map(sample, ((j, prime) for j in numbers), chunksize=2), 1
        ):
            a[i - 1] = row[cols] if saved else row
            if not saved and i in limits:
                with np.errstate(all="ignore"):
                    r, c, d = minor(a[:i], prime)
                assert len(r) <= 3017
                cert = dict(
                    prime=prime,
                    points=i,
                    rank_lower_bound=len(r),
                    rows=r,
                    columns=c,
                    determinant=d,
                )
                write(args.out / ("jacobi-%d.json" % i), cert)
                print(i, len(r), d, flush=True)
    if saved:
        with np.errstate(all="ignore"):
            d = square_determinant(a, prime)
        if prime == saved["prime"]:
            assert d == saved["determinant"]
        print(
            json.dumps(dict(prime=prime, size=len(numbers), determinant=d)), flush=True
        )


if __name__ == "__main__":
    main()
