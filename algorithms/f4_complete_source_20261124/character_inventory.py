"""Exact adjacent cochain dimensions for the F4 second incoming sources.

The colour and symmetric-group character routines are byte-for-byte public
inputs pinned in character_inputs.json.  This program does not construct
native source cochains or infer any differential rank from a character.
"""

import argparse
from fractions import Fraction
from functools import cache
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "character_backend"))
from b2_schur_complex import centralizer, character, multiplicity, partitions
from f4_singlet_characters import decomposition


def verify_inputs():
    records = json.loads((HERE / "character_inputs.json").read_text())
    for record in records:
        if sha256((HERE / record["path"]).read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError("character input digest mismatch: " + record["path"])
    return records


@cache
def singlet_trace(cycles):
    # Invariant inner product of two virtual Weyl characters.  Splitting the
    # Adams factors preserves the exact character while reducing memory.
    left, right = [], []
    for length in sorted(cycles, reverse=True):
        (left if sum(left) <= sum(right) else right).append(length)
    first = decomposition("F4", tuple(left))
    second = decomposition("F4", tuple(right))
    return sum(value * second.get(weight, 0) for weight, value in first.items())


def admitted_types(power, order):
    spin = (power - order, power - order)
    flavour = (3,) * order
    total = 2 * power + order
    signatures = [
        (
            tuple(sorted(first + second)),
            Fraction(
                character(spin, first) * character(flavour, second),
                centralizer(first) * centralizer(second),
            ),
        )
        for first in partitions(2 * (power - order))
        for second in partitions(3 * order)
    ]
    records = []
    for partition in partitions(total):
        shape = partition[::-1]
        coefficient = sum(
            value * character(shape, cycles) for cycles, value in signatures
        )
        if not coefficient:
            continue
        if coefficient.denominator != 1 or coefficient <= 0:
            raise ArithmeticError("nonintegral branching multiplicity")
        highest = (
            *shape[:2],
            *(max(0, sum(length >= column for length in shape) - 2)
              for column in (1, 2, 3)),
        )
        shifted = tuple(
            value + shift
            for value, shift in zip(highest, (order, order, -order, -order, -order))
        )
        dimensions = {
            str(length): multiplicity(shape, length, singlet_trace)
            for length in (total - 2, total - 1, total)
        }
        record = {
            "taylor_partition": shape,
            "branching_multiplicity": int(coefficient),
            "highest_weight_dotted_first": highest,
            "shifted_highest_weight": shifted,
            "cochain_dimensions_by_length": dimensions,
            "selected_by_full_taylor_ward_identity": shifted == (power, power, 0, 0, 0),
        }
        records.append(record)
        print(json.dumps({"progress": power, "order": order, **record}), flush=True)
    expected = (power - order, power - order) + (3,) * order
    selected = [r for r in records if r["selected_by_full_taylor_ward_identity"]]
    if len(selected) != 1 or tuple(selected[0]["taylor_partition"]) != expected:
        raise ArithmeticError("Taylor selection did not isolate the expected source")
    return {
        "power": power,
        "order": order,
        "native_charge_dotted_first": [power - order, power - order, order, order, order],
        "source_length": total - 1,
        "all_spin_flavour_singlet_taylor_types": records,
        "full_singlet_dimensions_by_length": {
            str(length): sum(
                record["branching_multiplicity"] * record["cochain_dimensions_by_length"][str(length)]
                for record in records
            )
            for length in (total - 2, total - 1, total)
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--power", type=int, choices=(6, 7), action="append")
    parser.add_argument("--order", type=int, choices=(1, 2), default=2)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = perf_counter()
    inputs = verify_inputs()
    records = [admitted_types(power, args.order) for power in args.power or (6, 7)]
    for record in records:
        if record["power"] == 6 and args.order == 2:
            if record["full_singlet_dimensions_by_length"]["13"] != 70238:
                raise ArithmeticError("sixth-power full-source control failed")
            selected = next(r for r in record["all_spin_flavour_singlet_taylor_types"]
                            if r["selected_by_full_taylor_ward_identity"])
            if list(selected["cochain_dimensions_by_length"].values()) != [10007, 3720, 704]:
                raise ArithmeticError("sixth-power selected-source control failed")
    output = {
        "schema": "f4-adjacent-character-inventory-v1",
        "arithmetic": "integer Weyl characters and rational symmetric-group characters",
        "scope": "Cochain dimensions only; native bases and differential ranks are not certified here.",
        "inputs": inputs,
        "records": records,
    }
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"output": str(args.out), "elapsed_seconds": perf_counter() - started}), flush=True)


if __name__ == "__main__":
    main()
