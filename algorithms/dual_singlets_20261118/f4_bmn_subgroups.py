"""Branch the two one-fermion rectangular F4 sources to reductive subalgebras."""
import argparse
from collections import Counter
from functools import lru_cache
import gzip
from itertools import permutations
import json
from pathlib import Path

from f4_bmn_root_levi import PRIME, verify
from su3_levi_branching import branches, sign
from su3_levi_kernel import check, key


UNIT = (0, 0, 0)


@lru_cache(None)
def weights(partition):
    first, second, third = partition
    multiplicities = Counter()
    for upper in range(second, first+1):
        for lower in range(third, second+1):
            for left in range(lower, upper+1):
                multiplicities[left, upper+lower-left, sum(partition)-upper-lower] += 1
    return multiplicities


@lru_cache(None)
def constituents(shape):
    """Pairs of GL(3) polynomials and their tensor multiplicity in shape."""
    parts = [(i, j, k) for i in range(shape[0]+1)
             for j in range(min(i, shape[1])+1) for k in range(min(j, shape[2])+1)]
    total = tuple(t+r for t, r in zip(shape, (2, 1, 0)))
    result = []
    for left in parts:
        shifted = tuple(t+r for t, r in zip(left, (2, 1, 0)))
        for right in parts:
            if sum(left)+sum(right) != sum(shape):
                continue
            mult = weights(right)
            copies = sum(sign(order)*mult[tuple(total[i]-shifted[order[i]] for i in range(3))]
                         for order in permutations(range(3)))
            assert copies >= 0
            if copies:
                result.append((left, right, copies))
    if shape[0] == shape[2]:
        assert all(copies == 1 and right == tuple(shape[0]-v for v in reversed(left))
                   for left, right, copies in result)
    return tuple(result)


def central_branches(width):
    outer = (width, width, width)
    rows = []
    for s in range(width+1):
        rows.append(((width, width, width-s), 1))
    for s in range(width):
        rows.append(((width, width-1, width-s-1), 0))
    for s in range(1, width+1):
        rows.append(((width-1, width-1, width-s), 0))
    general = [r for parity in (1, -1) for r in branches((*outer, 0, 0), 3*width-1, parity)]
    assert sorted((tuple(r["weight"][:3]), sum(r["weight"])-r["length"]) for r in general) == sorted(rows)
    assert all(r["copies"] == 1 for r in general)
    return rows


def evaluate(records, groups, width, centre):
    active = central_branches(width) if centre else [((width,)*3, 1)]
    missing = set()

    def simple(group, shape, defect):
        if shape == UNIT:
            return (0, 0, 0) if defect else (1, 1, 0)
        length = sum(shape)-defect
        if length <= 0:
            return 0, 0, 0
        weight = (*shape, 0, 0)
        label = key(weight, length)
        if label not in records[group]:
            missing.add((group, label))
            return 0, 0, 0
        return check(records[group][label], weight, length, PRIME, None if group == "SU2" else group)

    size = classical = kernel = 0
    for shape, defect in active:
        if len(groups) == 1:
            c, h, k = simple(groups[0], shape, defect)
            size += c
            classical += h
            kernel += k
            continue
        for left, right, copies in constituents(shape):
            for left_defect in range(defect+1):
                # The unit cannot carry the unique matter fermion.
                if (left == UNIT and left_defect) or (right == UNIT and defect-left_defect):
                    continue
                a, b, k = simple(groups[0], left, left_defect)
                x, y, l = simple(groups[1], right, defect-left_defect)
                size += copies*a*x
                classical += copies*b*y
                # ker(f tensor g) = ker(f) tensor V + U tensor ker(g).
                kernel += copies*(k*y+b*l)
    if missing:
        return dict(missing=[list(item) for item in sorted(missing)])
    return dict(cochains=size, classical_upper=classical, Cartan_kernel_upper=kernel)


def run(files):
    records = {}
    for group, path in files.items():
        data = json.loads(gzip.decompress(path.read_bytes()))
        widths = [4, 5] if key((5, 5, 5, 0, 0), 14) in data else [4]
        verify(data, group, widths)
        records[group] = data
    configurations = {"B3+T1": (("SO7",), True), "C3+T1": (("Sp3",), True),
                      "A2+A1+T1": (("SU3","SU2"), True),
                      "A2+A2": (("SU3","SU3"), False), "C3+A1": (("Sp3","SU2"), False)}
    return {str(6*width): {name:evaluate(records, groups, width, centre)
                          for name, (groups, centre) in configurations.items()}
            for width in (4, 5)}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for label in ("SU2", "SU3", "SO7", "Sp3"):
        p.add_argument("--"+label, required=True, type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("The character replay requires assertions")
    result = run({g:getattr(args, g) for g in ("SU2","SU3","SO7","Sp3")})
    args.out.write_text(json.dumps(result, indent=2)+"\n")
