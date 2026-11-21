"""Check the central-degree decomposition of the F4 level-30 BMN source on its rank-three Levis."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path

from f4_bmn_root_levi import PRIME, verify as verify_simple
from f4_bmn_subgroups import UNIT, central_branches, constituents
from su3_levi_branching import branches
from su3_levi_kernel import check, key
from trace_descent import dimensions


@lru_cache(None)
def remove_centres(weight, length, centres):
    if centres == 0:
        return {(weight[:3], length): 1}
    result = Counter()
    for parity in (-1, 1):
        for row in branches(weight, length, parity):
            for block, multiplicity in remove_centres(tuple(row["weight"]), row["length"], centres-1).items():
                result[block] += row["copies"] * multiplicity
    return result


def central_multiplicities(width, centres):
    """Cauchy for central scalars and vertical Pieri for the unique matter fermion."""
    def dimension(shape):
        if any(shape[centres:]):
            return 0
        parts = (tuple(shape) + (0,)*centres)[:centres]
        value = Fraction(1)
        for i, j in combinations(range(centres), 2):
            value *= Fraction(parts[i] - parts[j] + j - i, j - i)
        if value.denominator != 1:
            raise ArithmeticError("Nonintegral flavour multiplicity")
        return value.numerator

    answer = {}
    for active, central, copies in constituents((width,)*3):
        assert copies == 1
        size = sum(active)
        with_active_fermion = dimension(central)
        if size > 1 and with_active_fermion:
            answer[active, size-1] = with_active_fermion
        with_central_fermion = 0
        for i, j in combinations(range(3), 2):
            smaller = list(central)
            smaller[i] -= 1
            smaller[j] -= 1
            if smaller[0] >= smaller[1] >= smaller[2] >= 0:
                with_central_fermion += centres*dimension(smaller)
        if with_central_fermion:
            answer[active, size] = with_central_fermion
    return answer


def lower_rank(folder, weight, length):
    records = {}
    for group in ("SU2", "SU3", "Sp2"):
        records[group] = json.loads(gzip.decompress((folder / (group + ".json.gz")).read_bytes()))
        verify_simple(records[group], group, widths=(5,))

    @lru_cache(None)
    def simple(group, shape, defect):
        if shape == UNIT:
            return (1, 1, 0) if defect == 0 else (0, 0, 0)
        y = sum(shape) - defect
        if not y:
            return 0, 0, 0
        w = (*shape, 0, 0)
        return check(records[group][key(w, y)], w, y, PRIME, None if group == "SU2" else group)

    output = {}
    for name, groups, centres in (("B2+T2", ("Sp2",), 2), ("A2+T2", ("SU3",), 2),
                                  ("A1+A1+T2", ("SU2", "SU2"), 2), ("A1+T3", ("SU2",), 3)):
        dim = upper = kernel = 0
        multiplicities = central_multiplicities(weight[0], centres)
        if multiplicities != remove_centres(weight, length, centres):
            raise ValueError("Central Cauchy and Pieri multiplicities disagree with iterated branching")
        for (shape, y), mult in multiplicities.items():
            defect = sum(shape)-y
            if defect not in (0, 1):
                raise ValueError("Unexpected matter-fermion degree")
            if len(groups) == 1:
                c, h, k = simple(groups[0], shape, defect)
                dim += mult*c
                upper += mult*h
                kernel += mult*k
                continue
            for left, right, copies in constituents(shape):
                for i in range(defect + 1):
                    if (left == UNIT and i) or (right == UNIT and defect-i):
                        continue
                    c, h, k = simple(groups[0], left, i)
                    d, j, l = simple(groups[1], right, defect-i)
                    dim += mult*copies*c*d
                    upper += mult*copies*h*j
                    kernel += mult*copies*(k*j + h*l)
        output[name] = dict(cochains=dim, classical_upper=upper, cartan_kernel_upper=kernel)
    return output


def verify(folder):
    raw = gzip.decompress((folder / "level30_central_minors.json.gz").read_bytes())
    receipt = json.loads((folder / "level30_central_receipt.json").read_bytes())
    if hashlib.sha256(raw).hexdigest() != receipt["uncompressed_minors_sha256"]:
        raise ValueError("The central-sector minor archive differs from its receipt")
    records = json.loads(raw)
    total_length = 14
    weight = (5, 5, 5, 0, 0)
    even = branches(weight, total_length, 1)
    odd = branches(weight, total_length, -1)
    expected = sorted((tuple(row["weight"][:3]), sum(row["weight"]) - row["length"])
                      for row in even + odd)
    if expected != sorted(central_branches(5)) or any(row["copies"] != 1 for row in even + odd):
        raise ValueError("The abelian character and the rectangular branching disagree")
    active = [row for row in even if row["length"] < total_length]
    if [row["weight"] for row in even if row["length"] == total_length] != [weight]:
        raise ValueError("Unexpected pure simple summand")

    result = dict(physical_level=30, source_length=total_length, prime=PRIME,
                  odd_central_branches=len(odd), positive_even_central_branches=len(active),
                  pure_active_weight=weight, groups={})
    for group, size, colour in (("SO7", 7, "SO"), ("Sp3", 6, "Sp")):
        if set(records[group]) != {key(row["weight"], row["length"]) for row in active}:
            raise ValueError(f"Unexpected positive even central records for {group}")
        sectors = []
        totals = Counter()
        for row in sorted(active, key=lambda row: (-row["length"], row["weight"])):
            w, length = tuple(row["weight"]), row["length"]
            data = records[group][key(w, length)]
            dim, upper, kernel = check(data, w, length, PRIME, group)
            character = sum(dimensions((*w[3:], *w[:3]), length, size, colour).values())
            if dim != character:
                raise ValueError(f"Character mismatch for {group} at {w} and length {length}")
            counts = dict(cochains=dim, incoming=data["incoming"]["rank"],
                          outgoing=data["outgoing"]["rank"], stacked=data["outgoing_and_Cartan"]["rank"],
                          classical_upper=upper, cartan_kernel_upper=kernel)
            totals.update(counts)
            sectors.append(dict(central_letters=total_length-length, active_weight=w[:3],
                                active_length=length, **counts))
        result["groups"][group] = dict(sectors=sectors, totals=dict(totals))
    result["semisimple_rank_at_most_two"] = lower_rank(folder, weight, total_length)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("Enable assertions to replay the branching and minor checks")
    folder = Path(__file__).resolve().parent / "f4" / "levi"
    args.out.write_text(json.dumps(verify(folder), indent=2) + "\n")
