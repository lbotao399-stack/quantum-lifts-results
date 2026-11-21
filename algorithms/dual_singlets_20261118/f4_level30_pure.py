"""Check the two pure rank-three summands and restrict the level-30 F4 source to C3+A1."""
import argparse
import gzip
import json
from pathlib import Path

from f4_bmn_root_levi import PRIME
from f4_bmn_subgroups import UNIT, constituents
from su3_levi_kernel import check, key
from trace_descent import dimensions


FOLDER = Path(__file__).resolve().parent / "f4/levi"
WEIGHT = (5, 5, 5, 0, 0)


def read(name):
    return json.loads(gzip.decompress((FOLDER / name).read_bytes()))


def calculate():
    records = {g: read(g + "_level30_pure.json.gz") for g in ("SO7", "Sp3")}
    pure = {}
    for group, size, colour in (("SO7", 7, "SO"), ("Sp3", 6, "Sp")):
        row = records[group]
        cochains, upper, kernel = check(row, WEIGHT, 14, PRIME, group)
        sizes = [sum(dimensions((0, 0, 5, 5, 5), length, size, colour).values())
                 for length in (13, 14, 15)]
        if sizes != [row["parent"], cochains, row["next"]]:
            raise ValueError(f"Independent colour character disagrees for {group}")
        incoming, outgoing, stacked = (row[k]["rank"] for k in ("incoming", "outgoing", "outgoing_and_Cartan"))
        if outgoing != sizes[-1] or stacked - outgoing != upper or kernel != 0:
            raise ValueError(f"The bounds on the rational cohomology do not coincide for {group}")
        pure[group] = dict(adjacent_cochains=sizes, incoming=incoming, outgoing=outgoing,
                           outgoing_and_Cartan=stacked, classical_dimension=upper, Cartan_kernel=kernel)

    rank_one = read("SU2.json.gz")
    symplectic = read("Sp3.json.gz")
    symplectic.update(read("level30_central_minors.json.gz")["Sp3"])
    symplectic[key(WEIGHT, 14)] = records["Sp3"]
    branches = []
    zero = 0
    for left, right, copies in constituents(WEIGHT[:3]):
        for left_defect in (0, 1):
            right_defect = 1 - left_defect
            if (left == UNIT and left_defect) or (right == UNIT and right_defect):
                continue
            if right == UNIT:
                right_upper = 1
            else:
                length = sum(right) - right_defect
                if length <= 0:
                    zero += copies
                    continue
                weight = (*right, 0, 0)
                _, right_upper, right_kernel = check(rank_one[key(weight, length)], weight, length, PRIME)
                if right_kernel:
                    raise ValueError("A rank-one factor still has a Cartan kernel")
            if not right_upper:
                zero += copies
                continue
            length = sum(left) - left_defect
            if left == UNIT:
                upper, kernel = 1, 0
            elif length <= 0:
                upper, kernel = 0, 0
            else:
                weight = (*left, 0, 0)
                _, upper, kernel = check(symplectic[key(weight, length)], weight, length, PRIME, "Sp3")
            branches.append(dict(symplectic_weight=left, symplectic_length=length,
                                 symplectic_upper=upper, symplectic_kernel=kernel,
                                 rank_one_weight=right, rank_one_length=sum(right)-right_defect,
                                 rank_one_upper=right_upper, multiplicity=copies))
    return dict(level=30, length=14, prime=PRIME, pure=pure,
                C3_A1=dict(branches_with_zero_rank_one_cohomology=zero, other_branches=branches,
                           classical_upper=sum(b["multiplicity"]*b["rank_one_upper"]*b["symplectic_upper"] for b in branches),
                           Cartan_kernel_upper=sum(b["multiplicity"]*b["rank_one_upper"]*b["symplectic_kernel"] for b in branches)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Enable assertions for the character checks")
    args.out.write_text(json.dumps(calculate(), indent=2) + "\n")
