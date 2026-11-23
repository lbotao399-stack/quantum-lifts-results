"""Independently audit four native polynomial Bianchi systems over the rationals.

The exterior Lie bracket and polynomial matrix assembly below do not use the
baseline's raw_bracket. Lie structure constants, leading frames, and scalar
curves are shared baseline inputs. Ranks use python-flint in both implementations.
This verifies the local mark calculation, not a global cohomology or quantum lift.

Usage:
    python independent_native_verify.py /path/to/pinned/first_page

By default the deterministic result is compared with independent_audit.json.
Use --out PATH to write a newly generated receipt instead.
"""

import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

from flint import fmpq_mat


BASELINE_REPOSITORY = "https://github.com/lbotao399-stack/quantum-lifts-results"
BASELINE_COMMIT = "0b2e1a6724361fb2666089bfbe0c1a36a23d48a1"
BASELINE_DIRECTORY = "algorithms/dual_singlets_20261118/f4/j36_source/first_page"
BASELINE_SHA256 = {
    "nilpotent_tangent.py": "c0357a728c9087286ba2990421d1c08f5e9a68fa85ff8ae9824b940fde4943f8",
    "nilpotent_quadratic.py": "0d146798d0011a2028cd07a1ecc6693c3683836e806d72ffb8289a1025527040",
    "nilpotent_cubic.py": "5b3d8b65ee1f08bf38c2bd85da3767aa2c9a88afae7248b9d1d1e8340c05619d",
    "nilpotent_arcs.py": "7fc7dc2425b495c7780482e61ffe2962f4e6c1694dc6f93583b590727bf7e9c2",
    "nilpotent_curved_arcs.py": "9747a9c3693e3b8906f798b663d89e5fcea944a280222c363c0a5809624b0e37",
    "nilpotent_source_bounds.json.gz": "0b21a36647b31681c995eb34d9f0ba4368fb96d474c42124b2bbdb522eece087",
}


def independent_bracket(colour, left, right):
    """Exterior product with native Lie bracket, using Python Fraction only."""
    result = defaultdict(Fraction)
    for (axes, first), left_value in left.items():
        for (other, second), right_value in right.items():
            if set(axes).intersection(other):
                continue
            inversions = sum(x > y for x in axes for y in other)
            coefficient = (
                Fraction(str(left_value))
                * Fraction(str(right_value))
                * (-1 if inversions % 2 else 1)
            )
            for target, value in colour.br[first][second].items():
                result[(tuple(sorted(axes + other)), target)] += (
                    coefficient * Fraction(str(value))
                )
    return {key: value for key, value in result.items() if value}


def independent_polynomial_bracket(colour, left, right):
    result = defaultdict(Fraction)
    for first_order, first in enumerate(left):
        for second_order, second in enumerate(right):
            for key, value in independent_bracket(colour, first, second).items():
                result[(first_order + second_order, key)] += value
    return {key: value for key, value in result.items() if value}


def audit(baseline):
    if not __debug__:
        raise RuntimeError("Run without -O: verification assertions must be enabled.")
    baseline = baseline.resolve()
    for name, expected in BASELINE_SHA256.items():
        actual = sha256((baseline / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"Pinned baseline SHA-256 mismatch: {name}")
    sys.path.insert(0, str(baseline))
    from nilpotent_quadratic import calculate
    from nilpotent_arcs import linear_combination
    from nilpotent_curved_arcs import Differential, curve_for
    from nilpotent_tangent import simple_coordinates
    from native_bianchi import simultaneous_marks

    with gzip.open(baseline / "nilpotent_source_bounds.json.gz") as source:
        archive = json.load(source)
    records = []
    for quartet in ((8, 19, 25, 46), (8, 19, 27, 47)):
        _, frames, _, _, _ = calculate(quartet, 1, 7, True)
        old = next(
            row for row in archive["rational_tests"]
            if tuple(row["root_colours"]) == quartet
            and row["short_axis"] == 1
            and row["scalar_degree"] == 7
        )
        differential = Differential(
            frames.colour, simple_coordinates(frames.colour), quartet
        )
        seen = set()
        tested = 0
        candidates = (
            old["curved"]["additional_evaluations"]
            + old["straight"]["independent_evaluations"]
        )
        for record in candidates:
            tangent = tuple(tuple(item) for item in record["tangent"])
            if tangent in seen:
                continue
            seen.add(tangent)
            curve = curve_for(
                frames.colour,
                differential,
                linear_combination(frames.forms[1], dict(tangent)),
                3,
            )
            assert curve is not None
            scalar = [{((i,), c): 1 for i, c in enumerate(quartet)}] + curve[1:]
            assert not independent_polynomial_bracket(frames.colour, scalar, scalar)
            allowed, certificates, receipt = simultaneous_marks(
                frames.colour, quartet, curve, frames.forms[2], 3
            )
            for series in certificates:
                assert not independent_polynomial_bracket(frames.colour, scalar, series)

            # Rebuild every native column independently. There are 936 positive-
            # order unknowns; no coefficient equation is sampled or discarded.
            columns = []
            for order in range(1, 4):
                for axes in combinations(range(4), 2):
                    for colour in range(52):
                        series = [{} for _ in range(4)]
                        series[order] = {(axes, colour): 1}
                        columns.append(
                            independent_polynomial_bracket(frames.colour, scalar, series)
                        )
            assert len(columns) == 936
            for _, form in frames.forms[2]:
                columns.append(
                    independent_polynomial_bracket(frames.colour, scalar, [form])
                )
            rows = sorted(set().union(*(column.keys() for column in columns)))
            matrix = fmpq_mat([
                [str(column.get(row, 0)) for column in columns] for row in rows
            ])
            higher = fmpq_mat([
                [matrix[i, j] for j in range(936)] for i in range(matrix.nrows())
            ])
            full_rank = matrix.rank()
            higher_rank = higher.rank()
            leading_dimension = len(frames.forms[2]) - full_rank + higher_rank
            assert full_rank == receipt["rank"]
            assert leading_dimension == allowed.ncols()
            assert len(rows) == receipt["native_equation_rows"]
            records.append({
                "root_colours": list(quartet),
                "curve_index": tested,
                "native_equation_rows": len(rows),
                "leading_frame_dimension": len(frames.forms[2]),
                "higher_unknowns": 936,
                "independent_full_rank": full_rank,
                "independent_higher_rank": higher_rank,
                "allowed_leading_dimension": leading_dimension,
                "independent_fraction_mark_checks": len(certificates),
                "all_scalar_coefficients_zero": True,
                "all_mark_coefficients_zero": True,
            })
            tested += 1
            if tested == 2:
                break
        assert tested == 2
    return {
        "baseline_repository": BASELINE_REPOSITORY,
        "baseline_commit": BASELINE_COMMIT,
        "baseline_directory": BASELINE_DIRECTORY,
        "baseline_sha256": BASELINE_SHA256,
        "scope": "Four fixed native scalar curves; supplied leading frames; mark degree at most three. No global cohomology or quantum-lift claim.",
        "shared_inputs": "Pinned Lie structure constants, leading frames and scalar curves; python-flint rank backend.",
        "independent_operations": "Python Fraction exterior Lie bracket, every polynomial coefficient, complete native matrix assembly, full-versus-higher rank count.",
        "curves": records,
        "total_scalar_curves_verified": len(records),
        "total_mark_completions_verified": sum(
            row["independent_fraction_mark_checks"] for row in records
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Run without -O: verification assertions must be enabled.")
    result = audit(args.baseline)
    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n")
    else:
        expected = json.loads(Path(__file__).with_name("independent_audit.json").read_text())
        assert result == expected
    print(json.dumps({
        "scalar_curves_verified": result["total_scalar_curves_verified"],
        "mark_completions_verified": result["total_mark_completions_verified"],
        "allowed_leading_dimensions": [
            row["allowed_leading_dimension"] for row in result["curves"]
        ],
        "status": "PASS",
    }))


if __name__ == "__main__":
    main()
