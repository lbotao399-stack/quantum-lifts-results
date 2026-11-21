"""Compute the Taylor cochains available to the G2 level-36 BMN primary."""
import argparse
import json
from pathlib import Path

from b2_schur_complex import multiplicity, partner, transpose
from bmn_primary_index import g2_trace


def calculate(weight, order):
    shape = partner(weight, order)
    if shape is None:
        return dict(order=order, shape=None, cochains_by_length={}, euler=0)
    if order == 0 and shape != transpose(weight):
        raise ArithmeticError("The BMN source diagram does not reproduce its flavour partition")
    dimensions = {length: multiplicity(shape, length, g2_trace)
                  for length in range(1, sum(shape)+1)}
    return dict(order=order, shape=shape, cochains_by_length=dimensions,
                euler=sum((-1)**length * count for length, count in dimensions.items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, choices=range(6), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = calculate((6, 6, 6), args.order)
    args.out.write_text(json.dumps(result, indent=2)+"\n")
