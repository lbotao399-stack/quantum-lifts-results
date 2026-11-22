"""Verify the two rational second transposes in the SO(14), SO(15) and SO(16) quadratic targets."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from functools import cache
import gzip
import json
from pathlib import Path

from flint import fmpq_mat
from replay import pivots
from auxiliary import quotient
from relations import determinant_evaluation
from verify import encode, first_duals, path_types, pfaffian, pullback, unit, words

HERE = Path(__file__).resolve().parent


@cache
def inputs(colours, power):
    length = 2 * power + 2
    assert not path_types(power, power, 2, length - colours - 1)
    first = path_types(power - 1, power - 1, 5, length - colours)
    ((parts, plus, minus, _),) = first
    assert all(v == (False, "mp") for v in parts) and plus == minus == 3
    spec = dict(
        paths=["mp"] * len(parts),
        matrices=[unit(4, 3)] * len(parts),
        scalar_directions=[[]] * len(parts),
        single=(tuple(range(3)),) * 2,
    )
    forms, _ = first_duals(colours, (power, power, 2, 0, 0), range(3))
    first_rows = []
    for axes in ((0, 0), (1, 2)):
        rel = pfaffian(spec, axes, 1, length + 1, colours)
        first_rows.append(
            [sum(c * rel.coefficient(encode(w)) for w, c in f.items()) for f in forms]
        )
        rel.clear()
    assert first_rows[0] == [2 * x for x in first_rows[1]]
    a, b, c = first_rows[0]
    assert a and b == -2 * a and c == -(2 * power + 4) * a
    duals, checks = [], []
    for covector in ((2, 1, 0), (2 * power + 4, 0, 1)):
        f = Counter()
        for s, form in zip(covector, forms):
            f.update({w: s * c for w, c in form.items()})
        second, info = pullback(colours, {w: c for w, c in f.items() if c}, -1)
        duals.append([(encode(z), c) for z, c in second.items()])
        checks.append(info["checked_columns"])
    indices = json.loads((HERE / "result.json").read_text())["site_indices"]
    basis = words((5, 5, 4, 2, 2))
    extra = ((1, 1, 0, 0, 0),) * (power - 7)
    sites = [
        tuple(sorted(basis[j] + extra, key=lambda x: (sum(x), x))) for j in indices
    ]
    assert len(set(sites)) == 20 and all(len(w) > colours // 2 for w in sites)
    shapes = path_types(power - 2, power - 2, 8, length + 1 - colours)
    tasks = [(i, j) for i in range(len(shapes)) for j in quotient(i, shapes)[0]]
    assert len(shapes) == 46 and len(tasks) == 208
    assert sum(len(quotient(i, shapes)[1]) for i in range(len(shapes))) == 1848
    return length, shapes, tasks, duals, sites, first_rows, checks


def row(specification):
    colours, power, i, j, native = specification
    length, shapes, _, duals, sites, _, _ = inputs(colours, power)
    ordered, plane, rel = determinant_evaluation(
        i,
        j,
        native,
        shapes=shapes,
        length=length + 2,
        colours=colours,
        compress_pairs=not native,
    )
    scores = [str(sum(c * ordered(w) for w, c in functional)) for functional in duals]
    columns = [str(plane(w)) for w in sites]
    rel.clear()
    return [i, j, scores, columns]


def verify(colours, power, rebuild, workers, native):
    length, shapes, tasks, duals, sites, first, checks = inputs(colours, power)
    rows = json.loads(
        gzip.decompress((HERE / f"rows-{colours}-{power}.json.gz").read_bytes())
    )
    assert sorted((i, j) for i, j, _, _ in rows) == sorted(tasks)
    matrix = fmpq_mat([coord for _, _, _, coord in rows])
    values = fmpq_mat([scores for _, _, scores, _ in rows])
    piv = pivots(matrix.transpose())
    assert len(piv) == 20
    minor = fmpq_mat([rows[i][3] for i in piv])
    solution = minor.solve(fmpq_mat([rows[i][2] for i in piv]))
    assert minor.det() and matrix * solution == values
    selected = range(len(rows)) if rebuild else [piv[0], piv[1]]
    jobs = [(colours, power, *rows[i][:2], native) for i in selected]
    if workers == 1:
        for i, task in zip(selected, jobs, strict=True):
            assert row(task) == rows[i], rows[i][:2]
    else:
        with ProcessPoolExecutor(workers) as pool:
            for i, r in zip(selected, pool.map(row, jobs), strict=True):
                assert r == rows[i], rows[i][:2]
    return dict(
        colours=colours,
        power=power,
        first=[list(map(str, r)) for r in first],
        second_terms=list(map(len, duals)),
        second_columns=checks,
        path_types=len(shapes),
        zero_weight=1848,
        independent=len(rows),
        plane_rank=len(piv),
        augmented_rank=len(piv),
        pivot_rows=piv,
        minor=str(minor.det()),
        descent_coefficients=[[str(x) for x in s] for s in solution.tolist()],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--native", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--colours", type=int, choices=(14, 15, 16), nargs="+", default=(14, 15, 16)
    )
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("Enable assertions and use a positive worker count")
    print(
        json.dumps(
            [
                verify(n, n - 6, args.rebuild, args.workers, args.native)
                for n in args.colours
            ],
            indent=2,
        )
    )
