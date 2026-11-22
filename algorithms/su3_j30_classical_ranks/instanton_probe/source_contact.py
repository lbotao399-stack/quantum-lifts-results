"""Test source-current flavour selection and separate the three tree contact variations."""

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
import json
from random import Random
from cubic_colours import source
from source_fluctuations import add, trace

# Each coefficient carries six exterior test modes and at most one of three
# independent odd parameters for the resolved supersymmetry variation.
MARKERS = 7 << 6


@cache
def product_mask(a, b):
    if a & b or a & MARKERS and b & MARKERS:
        return 0, 0
    inv = 0
    x = a
    while x:
        bit = x & -x
        inv += (b & (bit - 1)).bit_count()
        x ^= bit
    return a | b, (-1) ** inv


def wedge(a, b):
    out = Counter()
    for x, v in a.items():
        for y, w in b.items():
            mask, sign = product_mask(x, y)
            if sign:
                out[mask] += sign * v * w
    return {m: v for m, v in out.items() if v}


def multiply(x, y):
    return [
        [add(*(wedge(x[i][k], y[k][j]) for k in range(3))) for j in range(3)]
        for i in range(3)
    ]


def traceless(rng):
    mat = [[rng.randint(-2, 2) for _ in range(3)] for _ in range(3)]
    mat[2][2] = -mat[0][0] - mat[1][1]
    return mat


def resolved_variations(seed):
    rng = Random(seed)
    fields = [
        [[{0: x} if x else {} for x in row] for row in traceless(rng)] for _ in range(3)
    ]
    for flavour in range(3):
        matrices = [traceless(rng) for _ in range(6)]
        left, right = fields[(flavour + 1) % 3], fields[(flavour + 2) % 3]
        ab, ba = multiply(left, right), multiply(right, left)
        matrix = []
        for i in range(3):
            row = []
            for j in range(3):
                entry = {1 << k: m[i][j] for k, m in enumerate(matrices) if m[i][j]}
                commutator = ab[i][j].get(0, 0) - ba[i][j].get(0, 0)
                if commutator:
                    entry[1 << (6 + flavour)] = commutator
                row.append(entry)
            matrix.append(row)
        assert not trace(matrix)
        fields.append(matrix)

    @cache
    def word(letters):
        if len(letters) == 1:
            return fields[letters[0]]
        return multiply(word(letters[:-1]), fields[letters[-1]])

    @cache
    def traced(letters):
        return trace(word(letters))

    combined = Counter()
    for coefficient, traces in source():
        term = {0: Fraction(coefficient)}
        for t in traces:
            term = wedge(term, traced(tuple(t)))
            if not term:
                break
        combined.update(term)
    rows = {
        m: [combined.get(m | (1 << (6 + i)), 0) for i in range(3)]
        for m in range(64)
        if m.bit_count() == 4
    }
    assert all(sum(v) == 0 for v in rows.values())
    for i, j in combinations(sorted(rows), 2):
        v, w = rows[i], rows[j]
        minor = v[0] * w[1] - v[1] * w[0]
        if minor:
            return dict(
                seed=seed,
                masks=[hex(i), hex(j)],
                coordinates=[list(map(str, v)), list(map(str, w))],
                minor=str(minor),
            )
    raise AssertionError(
        "The evaluation does not separate the two independent contact terms"
    )


def verify():
    # Each chiral scalar is a pair of indices in the three matter flavours.
    matter = set(range(3))
    scalars = [matter - {i} for i in range(3)]
    # The bilinear Poincare current carries the conjugate fourth R index.
    # Its D(bar phi_4i) factors also have the two complementary matter indices.
    epsilon_values = [len(x | y) == 4 for x in scalars for y in scalars]
    fermion_values = [i == 3 for i in matter]
    assert not any(epsilon_values + fermion_values)
    return dict(
        scalar_propagators_zero=len(epsilon_values),
        mixed_chirality_flavour_deltas_zero=len(fermion_values),
        classical_contact_samples=[resolved_variations(seed) for seed in (7, 19)],
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
