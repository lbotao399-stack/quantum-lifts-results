"""Extend exact polynomial families at the two unresolved local F4 sites.

This is a local necessary-condition test, not a calculation of the global d2.
Every accepted family is checked in the original native colour equations.
"""

import argparse
import gzip
import json
from pathlib import Path
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("--site", choices=("25", "27"), required=True)
    parser.add_argument("--degree", type=int, default=5)
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.degree < 1 or args.samples < 1:
        parser.error("enable assertions and choose positive degree and sample count")
    sys.path.insert(0, str(args.baseline))
    from nilpotent_curved_arcs import calculate
    from nilpotent_source_bounds import invariant_bound

    archive = json.load(gzip.open(args.baseline / "nilpotent_source_bounds.json.gz"))
    site = (8, 19, 25, 46) if args.site == "25" else (8, 19, 27, 47)
    record = next(row for row in archive["rational_tests"]
                  if tuple(row["root_colours"]) == site
                  and row["short_axis"] == 1 and row["scalar_degree"] == 7)
    start = time.perf_counter()
    result = calculate(site, 1, 7, record["straight"], args.degree, args.samples)
    fixed = invariant_bound(site, 1, 7, record["straight"], result)
    data = {
        "baseline_commit": "0b2e1a6724361fb2666089bfbe0c1a36a23d48a1",
        "test": "native polynomial commuting and Bianchi equations",
        "previous_polynomial_bound": record["curved"]["global_source_symbol_upper_bound"],
        "previous_finite_invariant_bound": record["finite_invariants"]["bound"],
        "elapsed_seconds": time.perf_counter() - start,
        "result": result,
        "finite_invariants": fixed,
        "scope": "local degree-seven symbols only; no global quantum lift claim",
    }
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps({k: v for k, v in data.items() if k not in ("result", "finite_invariants")}, indent=2))
    print("new polynomial bound", result["global_source_symbol_upper_bound"])
    print("new finite invariant bound", fixed["bound"])


if __name__ == "__main__":
    main()
