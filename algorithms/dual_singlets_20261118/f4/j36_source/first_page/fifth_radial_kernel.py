"""Integrate the fourth flavour of the native F4 instanton quartic for arbitrary external colours."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations
import argparse
import gzip
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

from flint import fmpz_mat, nmod_mat

from fifth_instanton import calculate as endpoint


def sign(left, right):
    swaps = 0
    while left:
        bit = left & -left
        swaps += (right & (bit - 1)).bit_count()
        left ^= bit
    return (-1) ** swaps


def inputs():
    data = json.loads(Path(__file__).with_name("fifth_quartic.json").read_bytes())
    physical = endpoint()
    colours = physical["doublets"]
    pairs = [
        (colours.index(i), colours.index(j), 2 * b)
        for i, j, b in physical["alternating_form"]
    ]
    entries = {}
    for *slots, value in data["nonzero_native_components"]:
        for p in set(permutations(colours.index(i) for i in slots)):
            entries[p] = Fraction(value)
    assert len(entries) == 840
    return pairs, entries


def complements(pairs, active):
    choices = []
    for used in combinations(range(7), active // 2):
        omitted = set(used)
        mask = 0
        coefficient = 1
        for i, (a, b, w) in enumerate(pairs):
            if i not in omitted:
                pair = (1 << a) | (1 << b)
                coefficient *= w * sign(mask, pair)
                mask |= pair
        rest = ((1 << 14) - 1) ^ mask
        choices.append(
            (rest, coefficient * sign(mask, rest) * factorial(7 - active // 2))
        )
    return choices


def sextic(pairs, entries):
    by_fourth = defaultdict(list)
    for (*triple, i), v in entries.items():
        by_fourth[i].append((triple, v))
    values = Counter()
    for i, j, b in pairs:
        for a, x in by_fourth[i]:
            for d, y in by_fourth[j]:
                if any(r == s for r, s in zip(a, d)):
                    continue
                key = tuple((1 << r) | (1 << s) for r, s in zip(a, d))
                parity = (-1) ** sum(r > s for r, s in zip(a, d))
                v = 16 * x * y / b
                assert v.denominator == 1
                values[key] += parity * int(v)
    result = {k: v for k, v in values.items() if v}
    assert all(
        all(result.get(p) == v for p in permutations(k)) for k, v in result.items()
    )
    return result


def flattening(values):
    pairs = sorted({a for a, b, c in values})
    columns = sorted({(b, c) for a, b, c in values})
    data = [[values.get((a, b, c), 0) for b, c in columns] for a in pairs]
    prime = 1009
    reduced, rank = nmod_mat(data, prime).rref()
    picked = [
        next(j for j in range(len(columns)) if reduced[i, j]) for i in range(rank)
    ]
    assert rank == len(pairs) == 91
    det = fmpz_mat([[row[j] for j in picked] for row in data]).det()
    assert det
    return dict(
        rows=len(pairs),
        nonzero_columns=len(columns),
        rank=rank,
        selected_columns=[columns[i] for i in picked],
        integer_determinant=str(det),
    )


def test_coefficients(pairs, entries, values):
    by_one = defaultdict(list)
    for (i, j, k, l), x in entries.items():
        by_one[i].append((1 << j, 1 << k, 1 << l, x))
    seed = defaultdict(list)
    for key, v in values.items():
        for i in range(14):
            if key[0] & (1 << i):
                seed[i].append((key, v))
    pf = complements(pairs, 0)[0][1] // factorial(7)
    assert pf == -1024

    @cache
    def exterior(a, b, c, d):
        if not a:
            return Fraction(int(b == c == d == 0))
        bit = a & -a
        ans = 0
        for j, k, l, v in by_one[bit.bit_length() - 1]:
            if b & j and c & k and d & l:
                phase = sign(j, b ^ j) * sign(k, c ^ k) * sign(l, d ^ l)
                ans += phase * v * exterior(a ^ bit, b ^ j, c ^ k, d ^ l)
        return ans

    @cache
    def power(a, b, c):
        if not a:
            return Fraction(int(b == c == 0))
        bit = a & -a
        value = 0
        for (i, j, k), v in seed[bit.bit_length() - 1]:
            if a & i == i and b & j == j and c & k == k:
                phase = sign(i, a ^ i) * sign(j, b ^ j) * sign(k, c ^ k)
                value += phase * v * power(a ^ i, b ^ j, c ^ k)
        return value

    tests = []
    for r in (1, 2, 3):

        def mask(start):
            return sum((1 << i) | (1 << j) for i, j, _ in pairs[start : start + r])

        cases = [(mask(0),) * 3, (mask(7 - r),) * 3, (mask(0), mask(1), mask(2))]
        # A non-paired colour pattern also checks an off-diagonal tensor coordinate.
        for key in sorted(values):
            if len(set(key)) > 1:
                allowed = [
                    (1 << i) | (1 << j)
                    for i, j, _ in pairs
                    if not any(((1 << i) | (1 << j)) & c for c in key)
                ]
                if len(allowed) >= r - 1:
                    spectator = sum(allowed[: r - 1])
                    cases.append(tuple(c | spectator for c in key))
                    break
        for a, b, c in cases:
            direct = sum(x * exterior(a, b, c, d) for d, x in complements(pairs, 2 * r))
            reduced = Fraction(factorial(7 - r) * pf * power(a, b, c), 16**r)
            assert reduced == direct
            tests.append(
                dict(vertices=2 * r, masks=[a, b, c], coefficient=str(reduced))
            )
    assert all(
        any(x["coefficient"] != "0" and x["vertices"] == k for x in tests)
        for k in (2, 4, 6)
    )
    return tests


def calculate(coefficients_path=None):
    pairs, entries = inputs()
    values = sextic(pairs, entries)
    sparse = [[*key, v] for key, v in sorted(values.items())]
    raw = (json.dumps(sparse, separators=(",", ":")) + "\n").encode()
    if coefficients_path is None:
        saved = Path(__file__).with_name("fifth_radial_sextic.json.gz")
        assert raw == gzip.decompress(saved.read_bytes())
    else:
        coefficients_path.write_bytes(gzip.compress(raw, mtime=0))
    vals = list(values.values())
    return dict(
        convention="entries are 16 times the sextic in the divided-power formula",
        pairs=pairs,
        sparse_sha256=sha256(raw).hexdigest(),
        native_quartic_entries=len(entries),
        sextic_entries=len(vals),
        coefficient_histogram=sorted(Counter(vals).items()),
        flattening=flattening(values),
        direct_exterior_tests=test_coefficients(pairs, entries, values),
        four_vertex_radial_prefactor=-61440,
        six_vertex_triplet_radial_prefactor=-20480,
        matched_source_current="open",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument(
        "--coefficients",
        type=Path,
        help="write the sparse sextic instead of comparing the archived coefficients",
    )
    a = p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    data = calculate(a.coefficients)
    a.out.write_text(json.dumps(data, indent=2) + "\n")
    print(
        data["sextic_entries"],
        data["flattening"]["rank"],
        data["direct_exterior_tests"][-4:],
    )
