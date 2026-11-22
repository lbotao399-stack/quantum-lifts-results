"""Character and rank checks for the G2 BMN census and the level-24 product space."""
import argparse
from collections import Counter, defaultdict
from hashlib import sha1, sha256
from itertools import permutations
import json
from math import factorial, prod
from pathlib import Path
import sys

from b2_schur_complex import centralizer, partitions
from bmn_primary_index import cochains, flavour_characters


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def census():
    pinned = read("upstream_odd.json")
    assert pinned["upstream_commit"] == "4256f207986029540060fb44dd494746828f65ce"
    incoming, joined = Counter(), Counter()
    for name, expected_blob in pinned["git_blobs"].items():
        data = (HERE / name).read_bytes()
        assert sha1(f"blob {len(data)}\0".encode() + data).hexdigest() == expected_blob
        for line in data.splitlines():
            row = json.loads(line)
            assert row["group"] == "G2"
            if not row["check"]:
                continue
            weight = tuple(row["lam"])
            key = weight, sum(weight) - row["n"]
            incoming[key] = max(incoming[key], row["r_up"])
            joined[key] = max(joined[key], row["r_up"] + row["grav"])
    summary, blocks = defaultdict(Counter), {}
    for y, weight_spaces in cochains(14, "G2")[1].items():
        if y == 0:
            continue
        for weight, (dimension,) in weight_spaces.items():
            key = weight, y
            out = incoming[weight, y + 1]
            assert joined[key] + out == dimension, (key, dimension, joined[key], out)
            h = dimension - incoming[key] - out
            assert h >= 0 and (y % 2 == 0 or h == 0)
            summary[2 * sum(weight)].update(blocks=1, dimension=dimension, cohomology=h)
            blocks[key] = (dimension, incoming[key], out, h)
    assert len(blocks) == 617 and sum(v["cohomology"] for v in summary.values()) == 675
    assert [blocks[(6, 3, 3), y] for y in (9, 10, 11)] == [
        (733, 307, 426, 0), (629, 426, 201, 2), (226, 201, 25, 0)]
    raw = read("level26.json")
    covered = set()
    for record in raw.values():
        w = tuple(record["weight"])
        row, = record["blocks"]
        assert record["prime"] == 2147483647 and w[3:] == (0, 0) and row["length"] == 9
        assert blocks[w[:3], 9] == (row["dimension"], row["incoming"], row["outgoing"], 0)
        covered.add(w[:3])
    assert covered == {w for w, y in blocks if sum(w) == 13 and min(w) >= 2 and y == 9}
    assert len(covered) == 8
    return {j: dict(s) for j, s in sorted(summary.items())}


def cartan_dimension():
    # Traces of the six rotations in the reflection representation of W(G2).
    rotations = (2, 1, -1, -2, -1, 1)
    weight, rho, length = (6, 3, 3), (2, 1, 0), 10
    letters = flavour_characters(sum(weight))
    total = 0
    for cycle in partitions(length):
        summed = sum(prod(rotations[g * k % 6] for k in cycle) for g in range(6))
        summed += 6 * prod(2 if k % 2 == 0 else 0 for k in cycle)
        assert summed % 12 == 0
        highest = 0
        for perm in permutations(rho):
            shifted = tuple(n + a - b for n, a, b in zip(weight, rho, perm))
            sign = (-1)**sum(a < b for i, a in enumerate(perm) for b in perm[i+1:])
            highest += sign * letters(cycle).get(shifted, 0)
        total += (factorial(length) // centralizer(cycle)) * (summed // 12) * highest
    assert total == 0
    return total


def check_minors():
    records = read("minors.json")
    sizes = {
        "basis-2147483647": (629,426,428,[3,5]),
        "basis-1000033": (629,426,428,[3,5]),
        "P0-2147483647": (629,426,427,[0]),
        "P0-1000033": (629,426,427,[0]),
        "partner-t1": (3215,1331,1390,None),
        "partner-t3": (235,72,84,None),
    }
    assert records.keys() == sizes.keys()
    for label, (dimension, incoming, joined, cols) in sizes.items():
        record = records[label]
        assert (record["dimension"], record["boundary"], record["augmented"]) == (
            dimension, incoming, joined)
        if cols is not None:
            assert record["selected"] == cols
        assert len(set(record["selected"])) == joined-incoming
        for key, expected in (("boundary_minor", incoming), ("joined_minor", joined)):
            minor = record[key]
            assert len(set(minor["rows"])) == len(set(minor["columns"])) == expected
            assert min(minor["rows"]) >= 0 and max(minor["rows"]) < dimension
            assert min(minor["columns"]) >= 0
            assert max(minor["columns"]) < record["predecessor"] + (0 if key == "boundary_minor" else record["products"])
            assert 0 < minor["det"] < record["prime"]
    sources = read("partner-ranks.json")
    assert sources["0"] == [733, 307, 426] and sources["2"] == [2205, 923, 1282]
    assert sources["1"] == [3215,1331,1825] and sources["3"] == [235,72,151]
    assert all(sources[t][1]+sources[t][2] == sources[t][0] for t in ("0","2"))
    assert all(records["partner-t"+t]["augmented"]+sources[t][2] == sources[t][0] for t in ("1","3"))
    return {k: v["joined_minor"]["det"] for k,v in sorted(records.items())}


def verify_files():
    for line in (HERE / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        assert Path(name).name == name and sha256((HERE / name).read_bytes()).hexdigest() == digest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", type=Path, metavar="PINNED_REPO")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        raise ValueError("run with assertions and positive worker count")
    verify_files()
    answer = dict(census=census(), cartan_dimension=cartan_dimension(), joined_determinants=check_minors())
    if args.rebuild:
        from rebuild import regenerate
        actual = regenerate(args.rebuild.resolve(), args.workers)
        assert actual == read("minors.json")
    print(json.dumps(answer, indent=2))
