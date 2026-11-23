"""Run the five colour couplings on classical predecessors with projected subwords."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import json
from multiprocessing import get_context
from pathlib import Path
import numpy as np

from global_source_boundaries import HERE, HEIGHTS, colour, minor, write
from global_source_scan import square_determinant
from marked_jet_blocks import MarkEvaluation, old_row

# Complete columns in a cut require no alternation across its colour edge.
CUTS = (
    (4, 0, 3, 0),
    (4, 4, 0, 0),
    (4, 0, 3, 0),
    (0, 4, 3, 0),
    (4, 0, 2, 0),
    (4, 1, 3, 0),
    (4, 0, 3, 1),
    (4, 4, 0, 1),
    (3, 0, 3, 0),
    (4, 0, 3, 0),
    (4, 0, 0, 3),
    (4, 4, 0, 0),
    (4, 0, 3, 1),
    (4, 1, 0, 3),
    (4, 3, 0, 0),
    (1, 4, 3, 0),
)
WIDTH = 780
N = 44
BASE = HERE / "jacobi-3088.json"


def bank(rng, weights, count, marks):
    result = []
    for n in range(count):
        groups = [
            int(x)
            for x in rng.permutation(
                [i for i, k in enumerate(weights) for _ in range(k)]
            )
        ]
        if marks == 0:
            sizes = [1] * len(groups)
        elif marks == 1:
            sizes = [2] + [1] * (len(groups) - 2)
        else:
            sizes = (
                [3] + [1] * (len(groups) - 3)
                if n % 2 == 0
                else [2, 2] + [1] * (len(groups) - 4)
            )
        word = []
        start = 0
        for size in rng.permutation(sizes):
            word.append((("x", "y", "z")[size - 1], *groups[start : start + size]))
            start += size
        # Evaluating a projected proper subword creates a separate internal colour edge.
        for _ in range(n % 3 + 1):
            if len(word) < 3:
                break
            size = int(rng.integers(2, min(len(word) - 1, 5) + 1))
            offset = int(rng.integers(len(word) - size + 1))
            part = word[offset : offset + size]
            options = ["ad", "sym", "tr"]
            if size >= 3 or any(t[0] == "part" for t in part):
                options += ["skew"]
            if size >= 4 or any(t[0] == "part" for t in part):
                options += ["jo", "jo"]
            channel = options[int(rng.integers(len(options)))]
            word[offset : offset + size] = [("part", channel, tuple(part))]
        result.append(tuple(word))
    return tuple(result)


@cache
def recipe(case):
    rng = np.random.default_rng(4562311 + case)
    cuts = CUTS[case]
    complement = tuple(h - x for h, x in zip(HEIGHTS, cuts))
    # Each bank contains direct third-jet and two-mark precursors in
    # the left factor or one mark on each side.
    external = case % 2 == 0
    left = bank(rng, cuts, N, 1 if external else 2)
    right = bank(rng, complement, N, 1 if external else 0)
    return cuts, left, right, external


@cache
def schedules():
    return tuple(
        np.asarray(j, dtype=int)
        for j in json.loads((HERE / "colour_screen.json").read_bytes())["indices"]
    )


def point(item):
    number, p, selected = item
    with np.errstate(all="ignore"):
        e = MarkEvaluation(14337 + number, p)
        old = old_row(e, BASE)
        if selected is None:
            return np.concatenate(
                [
                    old,
                    *(
                        e.evaluate_bank(*recipe(c), indices)
                        for c, indices in enumerate(schedules())
                    ),
                ]
            )
        # The saved minor orders its columns by pivot sequence.
        selected = np.asarray(selected)
        out = np.empty(len(selected), dtype=np.uint16)
        first = np.flatnonzero(selected < 2748)
        out[first] = old[selected[first]]
        for c, indices in enumerate(schedules()):
            offset = 2748 + c * WIDTH
            take = np.flatnonzero((selected >= offset) & (selected < offset + WIDTH))
            if len(take):
                out[take] = e.evaluate_bank(
                    *recipe(c), indices[selected[take] - offset]
                )
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=("search", "replay", "controls", "screen"))
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--points", type=int, default=3060)
    ap.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    ap.add_argument("--certificate", type=Path)
    args = ap.parse_args()
    assert __debug__ and args.workers > 0
    p = args.prime
    if args.action == "screen":
        # Include only explicitly evaluated nonzero columns and shuffle before retaining a fixed subset.
        # Discarded columns are never used to bound the rank from above.
        with np.errstate(all="ignore"):
            e = MarkEvaluation(97115, 1009)
            choices = []
            for case in range(len(CUTS)):
                all_indices = np.arange(5 * N * N, dtype=int)
                values = e.evaluate_bank(*recipe(case), all_indices)
                live = np.flatnonzero(values)
                rng = np.random.default_rng(733620 + case)
                assert len(live) >= WIDTH, (case, len(live))
                choices.append(sorted(map(int, rng.choice(live, WIDTH, replace=False))))
            write(
                HERE / "colour_screen.json",
                dict(prime=1009, colour_seed=97115, width=WIDTH, indices=choices),
            )
        return
    if args.action == "controls":
        with np.errstate(all="ignore"):
            e = MarkEvaluation(14420, p, brackets=True)
            for case, indices in enumerate(schedules()):
                assert not e.evaluate_bank(*recipe(case), indices).any()
            x = np.random.default_rng(241).integers(0, p, (26, 26)).astype(float)
            pieces = [e.decorate(x, "tr") * pow(26, -1, p) % p] + [
                e.decorate(x, k) for k in ("ad", "jo", "sym", "skew")
            ]
            assert np.array_equal(sum(pieces) % p, x)
            for j, m in enumerate(pieces):
                projections = [e.decorate(m, "tr") * pow(26, -1, p) % p] + [
                    e.decorate(m, k) for k in ("ad", "jo", "sym", "skew")
                ]
                assert all(
                    np.array_equal(z, m) if i == j else not np.any(z)
                    for i, z in enumerate(projections)
                )
        print("All five projectors and all boundary split controls passed.", flush=True)
        return
    saved = (
        json.loads(args.certificate.read_bytes()) if args.action == "replay" else None
    )
    ids = saved["rows"] if saved else range(args.points)
    width = len(ids) if saved else 2748 + len(CUTS) * WIDTH
    rows = np.empty((len(ids), width), dtype=np.uint16)
    colour(p)
    limits = {min(len(ids), 2816), min(len(ids), 3018), len(ids)}
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        for n, row in enumerate(
            pool.map(
                point,
                ((i, p, saved["columns"] if saved else None) for i in ids),
                chunksize=2,
            ),
            1,
        ):
            rows[n - 1] = row
            if not saved and n in limits:
                with np.errstate(all="ignore"):
                    r, c, d = minor(rows[:n], p)
                assert len(r) <= 3017
                rec = dict(
                    points=n,
                    prime=p,
                    rank_lower_bound=len(r),
                    rows=r,
                    columns=c,
                    determinant=d,
                )
                write(Path(f"channels-{n}.json"), rec)
                print(n, len(r), d, flush=True)
    if saved:
        with np.errstate(all="ignore"):
            d = square_determinant(rows, p)
        if p == saved["prime"]:
            assert d == saved["determinant"]
        print(json.dumps(dict(size=len(ids), prime=p, determinant=d)), flush=True)


if __name__ == "__main__":
    main()
