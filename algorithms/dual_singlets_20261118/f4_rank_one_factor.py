"""Test restriction to a Cartan in the A1 factor of an F4 subgroup at levels 24 and 30."""
import argparse
import gzip
import json
from pathlib import Path

from f4_bmn_subgroups import UNIT, constituents
from su3_levi_kernel import check, compute, key
from trace_descent import dimensions


DATA = Path(__file__).resolve().parent / "f4/levi"


def branches(width):
    for left, right, copies in constituents((width,) * 3):
        for right_defect in (0, 1):
            if (left == UNIT and not right_defect) or (right == UNIT and right_defect):
                continue
            yield left, right, sum(right)-right_defect, copies


def evaluate(records, prime):
    levels = {}
    for width in (4, 5):
        rows = []
        for left, right, length, copies in branches(width):
            weight = (*right, 0, 0)
            if right == UNIT:
                size, classical, kernel = 1, 1, 0
            elif length == 0:
                size, classical, kernel = 0, 0, 0
            else:
                size, classical, kernel = check(records[key(weight, length)], weight, length, prime)
                if sum(dimensions((0, 0, *right), length, 2, "SU+").values()) != size:
                    raise ArithmeticError("Rank-one cochain dimension disagrees with its character")
            if kernel:
                raise ArithmeticError("Cartan restriction has an unresolved kernel")
            rows.append(dict(other_weight=left, rank_one_weight=right, rank_one_length=length,
                             multiplicity=copies, cochains=size, cohomology_upper=classical,
                             Cartan_kernel_upper=kernel))
        solo, = (r for r in rows if r["other_weight"] == UNIT)
        levels[str(6*width)] = dict(branches=len(rows), empty=sum(r["cohomology_upper"] == 0 for r in rows),
                                    sum_rank_one_dimensions=sum(r["cochains"] for r in rows),
                                    rank_one_alone=solo, Cartan_kernel_bound=0, summands=rows)
    return levels


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, help="Regenerate the rank-one minors at an independent prime")
    parser.add_argument("--work", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if not __debug__ or bool(args.engine) != bool(args.work):
        parser.error("Enable assertions and supply --work with --engine")
    original = json.loads(gzip.decompress((DATA/"SU2.json.gz").read_bytes()))
    result = evaluate(original, 2147483647)
    if args.engine:
        tasks = {(right, length) for width in (4, 5) for _, right, length, _ in branches(width) if length}
        cases = {"rectangles": {"active": [dict(weight=(*right, 0, 0), length=length)
                                           for right, length in sorted(tasks)]}}
        records = compute(cases, args.engine, args.work, 2, 1000033)
        if evaluate(records, 1000033) != result:
            raise ArithmeticError("Independent rank-one calculation disagrees")
    args.out.write_text(json.dumps(result, indent=2) + "\n")
