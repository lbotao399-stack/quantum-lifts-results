"""Branch a highest-weight cochain block to active SU(2) tensors and abelian spectators."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import json
from pathlib import Path


# (name, SU(3) flavour degrees, dotted-spin degrees, target length, colour parity).
CASES = (
    ("j32_first", (7, 3, 3, 1, 1), 11, 1),
    ("j34_first", (7, 4, 3, 1, 1), 11, 1),
    ("j34_second", (6, 3, 2, 2, 2), 10, 1),
    ("j36_third_from10", (5, 2, 2, 3, 3), 8, -1),
    ("j36_third_from11", (5, 2, 2, 3, 3), 9, 1),
    ("j36_third_from12", (5, 2, 2, 3, 3), 10, -1),
)


def sign(positions):
    inversions = sum(i > j for k, i in enumerate(positions) for j in positions[k + 1:])
    return (-1) ** inversions


def shifted_orbit(weight):
    for flavour in permutations(range(3)):
        for spin in ((3, 4), (4, 3)):
            order = flavour + spin
            shift = (2, 1, 0, 1, 0)
            yield sign(flavour) * sign(spin), tuple(weight[i] + shift[i] for i in order)


@lru_cache(None)
def spectators(cap, length):
    """Exact abelian jet monomials, including each exterior generator at most once."""
    counts = {(0, 0, 0, 0, 0, 0): 1}
    for flavour in product((0, 1), repeat=3):
        odd = sum(flavour) % 2 == 0
        for plus in range(cap[3] + 1):
            for minus in range(cap[4] + 1):
                step = flavour + (plus, minus, 1)
                if not any(step[:5]):
                    continue
                following = counts.copy()
                for degree, multiplicity in counts.items():
                    new = tuple(x + d for x, d in zip(degree, step))
                    while all(x <= c for x, c in zip(new, cap + (length,))):
                        following[new] = following.get(new, 0) + multiplicity
                        if odd:
                            break
                        new = tuple(x + d for x, d in zip(new, step))
                counts = following
    return counts


def branches(weight, length, parity):
    shifted = tuple(x + d for x, d in zip(weight, (2, 1, 0, 1, 0)))
    counts = spectators(shifted, length)
    rows = []
    for active in product(*(range(x + 1) for x in weight)):
        if not (active[0] >= active[1] >= active[2] and active[3] >= active[4]):
            continue
        multiplicities = Counter()
        for phase, translated in shifted_orbit(active):
            spectator = tuple(x - y for x, y in zip(shifted, translated))
            if min(spectator) < 0:
                continue
            for n in range(length + 1):
                if (-1) ** n == parity:
                    multiplicities[n] += phase * counts.get(spectator + (n,), 0)
        for n, copies in sorted(multiplicities.items()):
            assert copies >= 0
            active_length = length - n
            if not copies or active_length > sum(active) or (not active_length and any(active)):
                continue
            rows.append(dict(weight=active, length=active_length, copies=copies))
    return rows


def plan():
    output = {}
    for name, weight, length, parity in CASES:
        output[name] = dict(weight=weight, length=length, colour_parity=parity,
                            active=branches(weight, length, parity))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    records = plan()
    args.out.write_text(json.dumps(records, separators=(",", ":")) + "\n")
    for key, record in records.items():
        print(key, len(record["active"]))
