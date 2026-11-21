"""Extract a Taylor component from two even and two odd directions at finite colour."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

from b2_schur_complex import centralizer, character, jet_character, partitions, partner, transpose
from trace_descent import dimensions
from unitary_characters import adjoint_traces


def strips(spin, first_tail, second_tail):
    if spin < 2 or first_tail < second_tail or second_tail < 0:
        raise ValueError("invalid two-by-two hook")
    for first in range(spin - 1):
        for second in range(first + 1):
            even = (spin-2-second, spin-2-first)
            for k in range(min(first_tail-second_tail, first-second) + 1):
                odd = (first_tail+2+first-k, second_tail+2+second+k, 0)
                yield (-1)**(first+second), (*even, *odd)


def schur_colour_dimensions(shape, length, colours):
    """Independent character of invariant traceless adjoint tensors of fixed parity."""
    if length < 1 or length > sum(shape):
        raise ValueError("invalid letter number")
    totals = [Fraction(), Fraction()]
    for cycle in partitions(length):
        ordinary, conjugated = adjoint_traces(cycle, colours)
        taylor = sum(coefficient*character(shape, slots)
                     for slots, coefficient in jet_character(cycle, sum(shape)-length))
        signed = Fraction((-1)**(length-len(cycle)), 2*centralizer(cycle))*taylor
        totals[0] += signed*(ordinary+conjugated)
        totals[1] += signed*(ordinary-conjugated)
    if any(v.denominator != 1 or v < 0 for v in totals):
        raise ValueError("invalid cochain multiplicities from the symmetric group")
    return tuple(map(int, totals))


def calculate(weight, order, length, colours):
    shape = partner(tuple(weight), order)
    if not shape or len(shape) < 2 or shape[0] != shape[1] or max(shape[2:], default=0) > 2:
        raise ValueError("partner must occur on two even and two odd Taylor directions")
    tail = transpose(shape[2:]) + (0, 0)
    rows = []
    for index, (sign, degrees) in enumerate(strips(shape[0], tail[0], tail[1])):
        parts = [sum(dimensions(degrees, length, colours, group).values())
                 for group in ("SU+", "SU-")]
        rows.append(dict(index=index, sign=sign, even_degrees=degrees[:2], odd_degrees=degrees[2:],
                         even_cochains=parts[0], odd_cochains=parts[1]))
    counts = tuple(sum(row["sign"]*row[parity+"_cochains"] for row in rows)
                   for parity in ("even", "odd"))
    if counts != schur_colour_dimensions(shape, length, colours):
        raise ValueError("two-flavour subtraction disagrees with the independent Schur character")
    raw = (order, order, *(v-order for v in weight))
    full = [sum(dimensions(raw, length, colours, group).values()) for group in ("SU+", "SU-")]
    return dict(group=f"SU({colours})", source_weight=list(weight), order=order, target_length=length,
                partition=list(shape), full_weight_even_and_odd=full, schur_even_and_odd=counts, rows=rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--weight", nargs=3, type=int, default=(8,5,4))
    p.add_argument("--order", type=int, default=2)
    p.add_argument("--length", type=int, default=10)
    p.add_argument("--colours", type=int, default=3)
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()
    record = calculate(a.weight, a.order, a.length, a.colours)
    a.out.write_text(json.dumps(record, indent=2)+"\n")
    print(record["full_weight_even_and_odd"], record["schur_even_and_odd"])
