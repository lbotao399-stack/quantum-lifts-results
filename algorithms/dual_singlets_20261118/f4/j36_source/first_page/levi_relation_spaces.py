"""Exact finite-colour relation sizes in each ordinary rank-three Levi diagram."""

from collections import Counter
from functools import cache, partial
from itertools import product
import json
from math import comb, factorial
from pathlib import Path

from source_symmetry import character, multiplicity, partitions
from b2_schur_complex import transpose
from levi_source_minors import branches
from levi_polarizations import components


@cache
def fixed_matchings(cycles):
    total = 1
    for length, number in Counter(cycles).items():
        total *= sum(factorial(number) * length**k //
                     (factorial(number - 2*k) * 2**k * factorial(k))
                     for k in range(number // 2 + 1) if length % 2 == 0 or 2*k == number)
    return total


@cache
def pair_character(group, cycles):
    # In orthogonal Schur-Weyl duality the admissible shapes have even rows.
    # For Sp(6), transpose these shapes and include the sign character separately.
    if group == "stable":
        return fixed_matchings(cycles)
    degree = sum(cycles) // 2
    allowed = [tuple(2*k for k in reversed(part)) for part in partitions(degree)
               if (len(part) <= 7 if group == "SO7" else max(part, default=0) <= 3)]
    return sum(character(shape, cycles) for shape in allowed)


@cache
def colour_trace(group, cycles):
    groups = sorted(Counter(cycles).items())
    total = 0
    for twists in product(*(range(n + 1) for _, n in groups)):
        fundamental = []
        ways = 1
        for (length, number), odd in zip(groups, twists):
            fundamental.extend([length] * (2 * (number - odd)))
            fundamental.extend([2 * length] * odd)
            ways *= comb(number, odd)
        total += (-1)**sum(twists) * ways * pair_character(group, tuple(sorted(fundamental)))
    result, remainder = divmod(total, 2**len(cycles))
    assert remainder == 0
    return result


def calculate():
    rows = {}
    for group in ("stable", "SO7", "Sp3"):
        trace = partial(colour_trace, group)
        assert [trace((1,) * n) for n in (1, 2, 3)] == [0, 1, 1]
        assert trace((1, 2)) == -1
        entries = []
        for shape in branches(1):
            degree = sum(shape)
            dims = [multiplicity(transpose(shape), y, trace) for y in range(degree-2, degree+1)]
            entries.append(dict(shape=shape, parent_source_terminal=dims))
        rows[group] = entries
        if group == "stable":
            continue
        for stable, finite in zip(rows["stable"], entries):
            finite["finite_relations"] = [s - f for s, f in zip(stable["parent_source_terminal"],
                                                                  finite["parent_source_terminal"])]
            assert all(r >= 0 for r in finite["finite_relations"])
        # Reconstruct the complete physical dimensions of the two outstanding
        # diagrams' smallest polarizations from the individual Schur summands.
        cost = json.loads(Path(__file__).with_name("levi_polarizations.json").read_bytes())[group]
        by_weight = {tuple(r["weight"]): r["sizes"] for r in cost["dimensions"]}
        for index in (3, 8):
            weight = tuple(cost["diagrams"][index]["selected"])
            degree = sum(weight)
            reconstructed = [0, 0, 0]
            for partition in partitions(degree):
                shape = partition[::-1]
                if len(shape) > 3 and shape[3] > 2:
                    continue
                coefficient = components(shape).get(weight, 0)
                if coefficient:
                    for j, y in enumerate(range(degree - 2, degree + 1)):
                        reconstructed[j] += coefficient * multiplicity(transpose(shape), y, trace)
            assert reconstructed == by_weight[weight]
            entries[index]["independent_physical_control"] = dict(weight=weight, reconstructed=reconstructed)
    return rows


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(), indent=2))
