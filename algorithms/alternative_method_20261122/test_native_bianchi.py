"""Compare fixed-preimage and simultaneous native mark equations on F4 curves."""

import argparse
import gzip
import json
from math import prod
from pathlib import Path
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("--site", choices=("25", "27"), default="27")
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--limit", type=int, default=2)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.degree < 1 or args.limit < 1:
        parser.error("enable assertions and choose positive degree and curve limit")
    sys.path.insert(0, str(args.baseline))
    from nilpotent_arcs import linear_combination
    from nilpotent_curved_arcs import Differential, curve_for, mark_columns
    from nilpotent_quadratic import calculate, array
    from nilpotent_source_bounds import invariant_bound
    from nilpotent_tangent import simple_coordinates
    from flint import fmpq, fmpq_mat
    from native_bianchi import simultaneous_marks

    q = (8, 19, 25, 46) if args.site == "25" else (8, 19, 27, 47)
    archive = json.load(gzip.open(args.baseline / "nilpotent_source_bounds.json.gz"))
    old = next(r for r in archive["rational_tests"] if tuple(r["root_colours"]) == q
               and r["short_axis"] == 1 and r["scalar_degree"] == 7)
    _, frames, mons, incoming, outgoing = calculate(q, 1, 7, True)
    differential = Differential(frames.colour, simple_coordinates(frames.colour), q)
    candidates = old["curved"]["additional_evaluations"] + old["straight"]["independent_evaluations"]
    evaluation_records = list(old["curved"]["additional_evaluations"])
    native_new = []
    all_rows = outgoing.tolist()
    for row in old["straight"]["independent_evaluations"] + evaluation_records:
        values = dict(row["evaluation"])
        all_rows.append([fmpq(values.get(i, 0)) for i in range(len(mons))])
    known_rank = fmpq_mat(all_rows).rank()
    seen = set()
    records = []
    for entry in candidates:
        u = tuple(tuple(item) for item in entry["tangent"])
        if u in seen:
            continue
        seen.add(u)
        curve = curve_for(frames.colour, differential, linear_combination(frames.forms[1], dict(u)), args.degree)
        assert curve is not None
        start = time.perf_counter()
        fixed, _ = mark_columns(frames.colour, differential, curve, range(len(frames.forms[2])), args.degree, frames.forms[2])
        fixed_seconds = time.perf_counter() - start
        start = time.perf_counter()
        all_marks, certificates, receipt = simultaneous_marks(frames.colour, q, curve, frames.forms[2], args.degree)
        native_seconds = time.perf_counter() - start
        combined = array(all_marks.nrows(), [[m[i, j] for i in range(m.nrows())] for m in (all_marks, fixed) for j in range(m.ncols())])
        assert combined.rank() == all_marks.ncols()
        for j, mark_series in enumerate(certificates):
            values = [prod(dict(u).get(s, 0) for s in factors) * all_marks[m, j]
                      for factors, m in mons]
            assert not any((fmpq_mat([values]) * incoming).entries())
            new_rank = fmpq_mat(all_rows + [values]).rank()
            if new_rank == known_rank:
                continue
            assert new_rank == known_rank + 1
            known_rank = new_rank
            all_rows.append(values)
            entry = {
                "tangent": u,
                "mark": [(i, str(all_marks[i, j])) for i in range(all_marks.nrows()) if all_marks[i, j]],
                "scalar_degree": len(curve) - 1,
                "mark_degree": max(i for i, y in enumerate(mark_series) if y),
                "evaluation": [(i, str(value)) for i, value in enumerate(values) if value],
            }
            evaluation_records.append(entry)
            native_new.append({**entry,
                "scalar_coefficients": [[(axes, c, str(v)) for (axes, c), v in x.items()] for x in curve],
                "mark_coefficients": [[(axes, c, str(v)) for (axes, c), v in y.items()] for y in mark_series],
            })
        rec = {
            "tangent": u,
            "scalar_degree": len(curve) - 1,
            "mark_degree_bound": args.degree,
            "fixed_leading_dimension": fixed.ncols(),
            "native_leading_dimension": all_marks.ncols(),
            "fixed_seconds": fixed_seconds,
            "native_seconds": native_seconds,
            "receipt": receipt,
        }
        records.append(rec)
        print(json.dumps(rec), flush=True)
        if len(records) == args.limit:
            break
    bound = len(mons) - incoming.rank() - known_rank
    second = {"global_source_symbol_upper_bound": bound,
              "additional_evaluations": evaluation_records}
    fixed = invariant_bound(q, 1, 7, old["straight"], second)
    result = {"root_colours": q, "records": records,
              "previous_polynomial_bound": old["curved"]["global_source_symbol_upper_bound"],
              "previous_invariant_bound": old["finite_invariants"]["bound"],
              "polynomial_bound": bound,
              "finite_invariants": fixed,
              "new_native_evaluations": native_new}
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print("new constraints", len(native_new), "polynomial bound", bound,
          "finite invariant bound", fixed["bound"], flush=True)


if __name__ == "__main__":
    main()
