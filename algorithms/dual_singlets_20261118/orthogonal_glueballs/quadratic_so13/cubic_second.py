"""Rebuild the second dual on the three one-flavour cubic endpoints in 14 and 15 colours."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import gzip
import json
from pathlib import Path

from flint import fmpq_mat
from auxiliary import quotient
from relations import determinant_evaluation
from replay import pivots
from verify import encode, first_duals, path_types, pullback, target_coordinates, words

HERE = Path(__file__).resolve().parent
CASES = ((14, 5), (15, 7))


@cache
def data(colours, power):
    length = 2 * power + 4
    weight = (power + 1, power, 3, 0, 0)
    assert len(target_coordinates(weight)[0]) == 3
    assert not path_types(power + 1, power, 3, length - colours - 1)
    assert not path_types(power, power - 1, 6, length - colours)
    shapes = path_types(power - 1, power - 2, 9, length + 1 - colours)
    first_rows, first_checks = first_duals(colours, weight, range(3))
    second_rows, second_checks = [], []
    for functional in first_rows:
        row, checked = pullback(colours, functional, -1)
        second_rows.append([(encode(w), c) for w, c in row.items()])
        second_checks.append(checked["checked_columns"])
    basis = words((power - 1, power - 2, 5, 2, 2))
    return length, shapes, first_checks, second_rows, second_checks, basis


def rebuild(task):
    colours, power, k, j, sites = task
    length, shapes, _, duals, _, basis = data(colours, power)
    ordered, plane, relation = determinant_evaluation(
        k,
        j,
        native=True,
        shapes=shapes,
        length=length + 2,
        colours=colours,
        extra_axes=3,
    )
    try:
        scores = [str(sum(c * ordered(w) for w, c in f)) for f in duals]
        coordinates = [str(plane(basis[i])) for i in sites]
        return scores, coordinates
    finally:
        relation.clear()


def verify(colours, power, rows, workers, all_sites):
    _, shapes, first, duals, second, basis = data(colours, power)
    components = [quotient(k, shapes, 3) for k in range(len(shapes))]
    assert [(k, j) for k, (reps, _) in enumerate(components) for j in reps] == [
        (k, j) for k, j, _, _ in rows
    ]
    all_indices = [i for i, w in enumerate(basis) if len(w) > colours // 2]
    columns = sorted({int(i) for _, _, _, r in rows for i in r})
    assert set(columns) <= set(all_indices)
    plane = fmpq_mat([[r.get(str(i), "0") for i in columns] for _, _, _, r in rows])
    scores = fmpq_mat([v for _, _, v, _ in rows])
    sites = [columns[i] for i in pivots(plane)]
    pivot_rows = pivots(plane.transpose())
    selected = [[r.get(str(i), "0") for i in sites] for _, _, _, r in rows]
    minor = fmpq_mat([selected[i] for i in pivot_rows])
    fitted = minor.solve(fmpq_mat([rows[i][2] for i in pivot_rows]))
    assert minor.det() and fmpq_mat(selected) * fitted == scores
    check_sites = all_indices if all_sites else sites
    jobs = [(colours, power, k, j, check_sites) for k, j, _, _ in rows]
    with ProcessPoolExecutor(workers) as pool:
        for (actual, coords), (_, _, expected, r) in zip(
            pool.map(rebuild, jobs), rows, strict=True
        ):
            assert actual == expected and coords == [
                r.get(str(i), "0") for i in check_sites
            ]
    return dict(
        colours=colours,
        level=3 * (2 * power + 1) + 6,
        finite_shapes=len(shapes),
        zero_weight_components=sum(len(v) for _, v in components),
        independent_components=len(rows),
        first_native_columns=first,
        second_native_columns=second,
        second_dual_lengths=list(map(len, duals)),
        plane_rank=len(pivot_rows),
        augmented_rank=len(pivot_rows),
        nonzero_uncorrected_scores=[
            sum(scores[i, j] != 0 for i in range(len(rows))) for j in range(3)
        ],
        sites=sites,
        pivot_rows=pivot_rows,
        minor=str(minor.det()),
        fitted=[[str(v) for v in r] for r in fitted.tolist()],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--all-sites", action="store_true")
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("Enable assertions and choose a positive worker count")
    rows = json.loads(gzip.decompress((HERE / "cubic_second.json.gz").read_bytes()))
    print(
        json.dumps(
            [
                verify(n, p, rows[f"{n}-{p}"], args.workers, args.all_sites)
                for n, p in CASES
            ],
            indent=2,
        )
    )
