"""Interpolate only the complement of the previously certified 2,879 boundary columns."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import json
from multiprocessing import get_context
from pathlib import Path
import numpy as np

from colour_channel_search import point as projected_point
from global_source_boundaries import HEIGHTS, HERE, minor, write
from global_source_scan import square_determinant
from mixed_colour_edges import EdgeEvaluation, transition_colour

CUTS = (
    (4, 0, 3, 0),
    (4, 4, 0, 0),
    (3, 0, 3, 0),
    (4, 0, 2, 0),
    (4, 1, 3, 0),
    (4, 3, 0, 0),
    (4, 0, 2, 1),
    (4, 1, 2, 0),
    (3, 0, 3, 1),
    (3, 1, 3, 0),
    (4, 2, 0, 1),
    (4, 2, 1, 0),
    (4, 0, 3, 0),
    (4, 4, 0, 0),
    (3, 0, 3, 0),
    (4, 0, 2, 0),
)
WIDTH = 300
N = 28


def bank(rng, weights, n, marks):
    output = []
    for row in range(n):
        groups = [
            int(c)
            for c in rng.permutation(
                [j for j, w in enumerate(weights) for _ in range(w)]
            )
        ]
        if marks == 0:
            sizes = [1] * len(groups)
        elif marks == 1:
            sizes = [2] + [1] * (len(groups) - 2)
        elif row % 2:
            sizes = [3] + [1] * (len(groups) - 3)
        else:
            sizes = [2, 2] + [1] * (len(groups) - 4)
        parts = []
        start = 0
        for k in rng.permutation(sizes):
            parts.append((("x", "y", "z")[k - 1], *groups[start : start + k]))
            start += k
        crossing = int(rng.integers(3, min(6, len(parts))))
        tail = [
            (("adj" if j == 0 or rng.integers(2) else "def"), *t)
            for j, t in enumerate(parts[crossing:])
        ]
        output.append(tuple(parts[:crossing]) + ("cross",) + tuple(tail))
    return tuple(output)


@cache
def model(c):
    cuts = CUTS[c]
    other = tuple(h - i for h, i in zip(HEIGHTS, cuts))
    rng = np.random.default_rng(818176 + c)
    external = c % 3 != 1
    left = bank(rng, cuts, N, 1 if external else 2)
    right = bank(rng, other, N, 1 if external else 0)
    return cuts, left, right, external


@cache
def screens():
    return [
        np.asarray(v, dtype=int)
        for v in json.loads((HERE / "mixed_screen.json").read_bytes())["indices"]
    ]


def solve_left(base, test, p):
    """Solve W A = C mod p by panel elimination. Return W and the determinant of A."""
    n = len(base)
    a = np.concatenate((base.T, test.T), axis=1).astype(float)
    determinant = 1
    for start in range(0, n, 32):
        stop = min(n, start + 32)
        width = stop - start
        panel = a[start:, start:stop].astype(np.int64)
        factors = []
        for k in range(width):
            viable = np.flatnonzero(panel[k:, k])
            if not len(viable):
                raise ValueError("singular retained pivot matrix")
            i = k + int(viable[0])
            if i != k:
                a[[start + k, start + i]] = a[[start + i, start + k]]
                panel[[k, i]] = panel[[i, k]]
                for col in factors:
                    col[[k, i]] = col[[i, k]]
                determinant = -determinant
            lead = int(panel[k, k])
            determinant = determinant * lead % p
            col = panel[:, k] * pow(lead, -1, p) % p
            col[: k + 1] = 0
            panel = (panel - col[:, None] * panel[k]) % p
            factors.append(col)
        top = a[start:stop, stop:].copy()
        for k, col in enumerate(factors):
            top[k + 1 :] = (top[k + 1 :] - col[k + 1 : width, None] * top[k]) % p
        a[start:stop, start:stop] = panel[:width]
        a[start:stop, stop:] = top
        if stop < n:
            mult = np.stack([c[width:] for c in factors], axis=1).astype(float)
            a[stop:, stop:] = (a[stop:, stop:] - mult @ top) % p
    rhs = a[:, n:].copy()
    for start in reversed(range(0, n, 32)):
        stop = min(start + 32, n)
        for k in range(stop - 1, start - 1, -1):
            rhs[k] = (rhs[k] - a[k, k + 1 : stop] @ rhs[k + 1 : stop]) % p
            rhs[k] = rhs[k] * pow(int(a[k, k]), -1, p) % p
        rhs[:start] = (rhs[:start] - a[:start, start:stop] @ rhs[start:stop]) % p
    w = rhs.T
    assert np.array_equal(w @ base % p, test)
    return w, determinant


@cache
def previous():
    return json.loads((HERE / "channels-3060.json").read_bytes())


def base_point(task):
    i, p = task
    return projected_point((i, p, previous()["columns"]))


def new_point(task):
    i, p, picked = task
    with np.errstate(all="ignore"):
        e = EdgeEvaluation(14337 + i, p)
        if picked is None:
            return np.concatenate(
                [
                    e.block_edges(model(c), indices)
                    for c, indices in enumerate(screens())
                ]
            )
        picked = np.asarray(picked, dtype=int)
        out = np.empty(len(picked), dtype=np.uint16)
        for c, indices in enumerate(screens()):
            take = np.flatnonzero(picked // WIDTH == c)
            if len(take):
                out[take] = e.block_edges(model(c), indices[picked[take] % WIDTH])
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=("controls", "screen", "search", "replay"))
    ap.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--certificate", type=Path, default=HERE / "mixed-minor.json")
    args = ap.parse_args()
    assert __debug__ and args.workers > 0
    p = args.prime
    if args.action == "screen":
        e = EdgeEvaluation(36981, 1009)
        active = []
        for c in range(len(CUTS)):
            indices = np.arange(N * N)
            with np.errstate(all="ignore"):
                v = e.block_edges(model(c), indices)
            live = np.flatnonzero(v)
            rng = np.random.default_rng(9813124 + c)
            # All sampled columns are genuine derivatives. Pad with other distinct
            # columns if the screening point sees fewer than WIDTH.
            if len(live) >= WIDTH:
                choose = rng.choice(live, WIDTH, replace=False)
            else:
                choose = np.r_[
                    live,
                    rng.choice(
                        np.flatnonzero(v == 0), WIDTH - len(live), replace=False
                    ),
                ]
            active.append(sorted(map(int, choose)))
        write(HERE / "mixed_screen.json", dict(seed=36981, indices=active))
        return
    if args.action == "controls":
        with np.errstate(all="ignore"):
            e = EdgeEvaluation(16771, p, brackets=True)
            for c, idx in enumerate(screens()):
                assert not e.block_edges(model(c), idx).any()
            rng = np.random.default_rng(24682)
            for n in (11, 32, 37, 94):
                base = rng.integers(0, p, size=(n, n))
                base[:4, :4] = 0
                trial = rng.integers(0, p, size=(6, n))
                weights, d = solve_left(base, trial, p)
                assert d == square_determinant(base, p)
                new = rng.integers(0, p, size=(n + 6, 5))
                small = (new[n:] - weights @ new[:n]) % p
                r, c, det = minor(small, p)
                union = np.concatenate(
                    (
                        np.r_[base, trial[r]],
                        new[np.r_[np.arange(n), n + np.asarray(r)]][:, c],
                    ),
                    axis=1,
                )
                assert square_determinant(union, p) == d * det % p
        print("Mixed-edge Jacobi and Schur determinant controls passed.", flush=True)
        return
    old = previous()
    num = len(old["rows"])
    save = (
        json.loads(args.certificate.read_bytes()) if args.action == "replay" else None
    )
    tests = save["new_rows"] if save else list(range(5000, 5168))
    points = list(old["rows"]) + tests
    assert len(set(points)) == len(points) and num == 2879
    transition_colour(p)
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        base = np.array(
            list(pool.map(base_point, ((i, p) for i in points), chunksize=2))
        )
        with np.errstate(all="ignore"):
            weights, d = solve_left(base[:num], base[num:], p)
        assert d == ({1009: 431, 1013: 100}[p])
        print("retained", num, d, flush=True)
        picked = save["new_columns"] if save else None
        raw = np.asarray(
            list(pool.map(new_point, ((i, p, picked) for i in points), chunksize=2)),
            dtype=np.float64,
        )
    with np.errstate(all="ignore"):
        reduced = (raw[num:] - weights @ raw[:num]) % p
    if save:
        det = square_determinant(reduced, p)
        if p == save["prime"]:
            assert d * det % p == save["determinant"]
        print(
            json.dumps(
                dict(
                    prime=p,
                    size=len(points),
                    base_minor=d,
                    complement_minor=det,
                    determinant=int(d * det % p),
                )
            ),
            flush=True,
        )
    else:
        r, c, det = minor(reduced, p)
        assert num + len(r) <= 3017
        result = dict(
            prime=p,
            previous_rank=num,
            new_rows=[tests[k] for k in r],
            new_columns=c,
            complement_minor=det,
            rank=num + len(r),
            determinant=int(d * det % p),
        )
        write(args.certificate, result)
        print("extended rank", num + len(r), "determinant", d * det % p, flush=True)


if __name__ == "__main__":
    main()
