"""Verify the second-page endpoint on the flavour-highest orthogonal determinant relations."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gzip
import json
from flint import fmpq_mat
from auxiliary import HERE, SHAPES, quotient
from relations import determinant_evaluation
from second import make_second
from verify import encode, words


def pivots(matrix):
    reduced, rank = matrix.rref()
    return [
        next(j for j in range(reduced.ncols()) if reduced[i, j]) for i in range(rank)
    ]


def checked_row(task):
    k, j, dual, sites, native = task
    ordered, plane, relation = determinant_evaluation(k, j, native)
    score = sum(c * ordered(w) for w, c in dual)
    coords = {str(i): str(value) for i, w in sites if (value := plane(w))}
    relation.clear()
    return [k, j, str(score), coords]


def verify(workers, rebuild, native, components):
    rows = json.loads(gzip.decompress((HERE / "rows.json.gz").read_bytes()))
    raw_count = sum(len(quotient(k)[1]) for k in range(len(SHAPES)))
    inputs = [(k, j) for k in range(len(SHAPES)) for j in quotient(k)[0]]
    assert len(SHAPES) == 46 and raw_count == 1848 and len(inputs) == 208
    assert inputs == [(k, j) for k, j, _, _ in rows]
    basis = words((5, 5, 4, 2, 2))
    sites = [(i, w) for i, w in enumerate(basis) if len(w) > 6]
    assert len(basis) == 1327 and len(sites) == 953
    nonzero = sorted({int(i) for _, _, _, d in rows for i in d})
    assert set(nonzero) <= dict(sites).keys()
    matrix = fmpq_mat([[d.get(str(i), "0") for i in nonzero] for _, _, _, d in rows])
    values = fmpq_mat([[r[2]] for r in rows])
    cols = pivots(matrix)
    selected = pivots(matrix.transpose())
    minor = fmpq_mat([[matrix[i, j] for j in cols] for i in selected])
    weights = minor.solve(fmpq_mat([[values[i, 0]] for i in selected]))
    assert len(cols) == len(selected) == 20 and minor.det() != 0
    assert all(
        sum(matrix[i, j] * weights[a, 0] for a, j in enumerate(cols)) == values[i, 0]
        for i in range(len(rows))
    )
    dual = [(encode(w), Fraction(c)) for w, c in make_second().items()]
    sample = list(range(len(rows))) if rebuild else selected
    extra = []
    if components:
        for k in range(len(SHAPES)):
            representatives, raw = quotient(k)
            remaining = [i for i in range(len(raw)) if i not in representatives]
            extra.extend(
                (k, i) for i in (remaining if components == "all" else remaining[:1])
            )
    with ProcessPoolExecutor(max_workers=workers) as pool:
        tasks = [(*inputs[i], dual, sites, native) for i in sample]
        for i, computed in zip(sample, pool.map(checked_row, tasks), strict=True):
            assert computed == rows[i], inputs[i]
        small = [(nonzero[j], basis[nonzero[j]]) for j in cols]
        factors = [Fraction(str(v[0])) for v in weights.tolist()]
        tasks = ((k, j, dual, small, native) for k, j in extra)
        for k, j, score, plane in pool.map(checked_row, tasks):
            assert sum(
                v * Fraction(plane.get(str(i), 0)) for v, (i, _) in zip(factors, small)
            ) == Fraction(score), (k, j)
    print(
        json.dumps(
            dict(
                path_types=len(SHAPES),
                weight_zero=raw_count,
                representatives=len(rows),
                checked=len(sample),
                discarded_components_checked=len(extra),
                rank=len(cols),
                augmented_rank=len(cols),
                minor=str(minor.det()),
                ordered_second_terms=len(dual),
                site_indices=[nonzero[j] for j in cols],
                weights=[str(v[0]) for v in weights.tolist()],
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument(
        "--native",
        action="store_true",
        help="use the original trace logarithm with every cut",
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument(
        "--components",
        choices=("sample", "all"),
        help="also evaluate dependent auxiliary components",
    )
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("Enable assertions and request a positive worker count")
    verify(args.workers, args.rebuild, args.native, args.components)
