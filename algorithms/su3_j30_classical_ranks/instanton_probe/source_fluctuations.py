"""Expand the public SU(3) R-007 polynomial on a charge-one scalar fluctuation."""

from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import prod
from pathlib import Path
import gzip
import json


HERE = Path(__file__).resolve().parent
SOURCE = (
    HERE.parents[1]
    / "dual_singlets_20261118/taylor_hooks/su3_bmn_sources/dmk_o330_hw.json.gz"
)


def add(*terms):
    result = Counter()
    for term in terms:
        result.update(term)
    return {i: c for i, c in result.items() if c}


def scale(term, factor):
    return {i: c * factor for i, c in term.items() if c * factor}


@cache
def sign(left, right):
    inversions = 0
    while left:
        lowest = left & -left
        inversions += (right & (lowest - 1)).bit_count()
        left -= lowest
    return (-1) ** inversions


def wedge(left, right):
    result = Counter()
    for i, a in left.items():
        for j, b in right.items():
            if not i & j:
                result[i | j] += sign(i, j) * a * b
    return {i: c for i, c in result.items() if c}


def matrix_product(left, right):
    return [
        [add(*(wedge(left[i][k], right[k][j]) for k in range(3))) for j in range(3)]
        for i in range(3)
    ]


def trace(matrix):
    return add(*(matrix[j][j] for j in range(3)))


def profiles():
    # Mode order within each R flavour: lower eta_1, eta_2, nu, barred nu.
    def variable(flavour, slot, coefficient):
        return {1 << (4 * flavour + slot): coefficient}

    column = [
        [variable(a, 0, 4), variable(a, 1, 4), variable(a, 2, 1)] for a in range(3)
    ]
    row = [[variable(a, 1, -4), variable(a, 0, 4), variable(a, 3, 1)] for a in range(3)]
    fields = []
    for a, b in [(1, 2), (2, 0), (0, 1)]:
        scalar = [
            [
                add(
                    wedge(column[b][i], row[a][j]),
                    scale(wedge(column[a][i], row[b][j]), -1),
                )
                for j in range(3)
            ]
            for i in range(3)
        ]
        # The ADHM correction has support on the two instanton colours.
        # Trace zero fixes it relative to the dyad. Divide each scalar
        # by its transverse nu^B barred-nu^A component in the final ratio.
        correction = scale(trace(scalar), Fraction(-1, 2))
        for i in range(2):
            scalar[i][i] = add(scalar[i][i], correction)
        fields.append(scalar)
    for a in range(3):
        fields.append(
            [
                [
                    add(
                        column[a][i] if j == 0 else {},
                        scale(row[a][j], -1) if i == 1 else {},
                    )
                    for j in range(3)
                ]
                for i in range(3)
            ]
        )
    assert all(not trace(field) for field in fields)
    return fields


def response(terms, direction, *, split=False):
    fields = profiles()
    variation = [[{0: v} if v else {} for v in line] for line in direction]
    assert not trace(variation)
    if not split:
        for i in range(3):
            for j in range(3):
                fields[0][i][j] = add(fields[0][i][j], variation[i][j])
    fields.append(variation)
    identity = [[{0: 1} if i == j else {} for j in range(3)] for i in range(3)]

    @cache
    def ordered(word):
        if not word:
            return identity
        return matrix_product(ordered(word[:-1]), fields[word[-1]])

    @cache
    def trace_by_cuts(word):
        if not split:
            return {0: trace(ordered(word))}
        # Distribute the three external scalar slots explicitly.
        locations = [j for j, v in enumerate(word) if v == 0]
        polynomials = {}
        for mask in range(1 << len(locations)):
            count = mask.bit_count()
            if count > 3:
                continue
            replaced = list(word)
            for bit, j in enumerate(locations):
                if mask & (1 << bit):
                    replaced[j] = 6
            polynomials[count] = add(
                polynomials.get(count, {}), trace(ordered(tuple(replaced)))
            )
        return polynomials

    result = Counter()
    for coefficient, traces in terms:
        polynomial = {0: {0: Fraction(coefficient)}}
        for word in traces:
            following = {}
            for i, a in polynomial.items():
                for j, b in trace_by_cuts(tuple(word)).items():
                    if i + j <= 3 and (term := wedge(a, b)):
                        following[i + j] = add(following.get(i + j, {}), term)
            polynomial = following
            if not polynomial:
                break
        result.update(polynomial.get(3 if split else 0, {}))
    ordered.cache_clear()
    return {mask: value for mask, value in result.items() if value}


def radial_cubic():
    # Evaluate the six residual radial jets on the specified test matrices. The scalar colour
    # projector contracts an incoming entry (1,2) with entry (2,1).
    h = [[1, 0, 0], [0, 1, 0], [0, 0, -2]]
    b = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]
    a = [[h[i][j] + b[i][j] for j in range(3)] for i in range(3)]
    u = [[0, 0, 0], [1, 0, 0], [0, 0, 0]]
    matrices = dict(A=a, B=b, C=h, U=u)

    def tr(word):
        value = [[int(i == j) for j in range(3)] for i in range(3)]
        for letter in word:
            value = [
                [
                    sum(value[i][k] * matrices[letter][k][j] for k in range(3))
                    for j in range(3)
                ]
                for i in range(3)
            ]
        return sum(value[i][i] for i in range(3))

    traces = json.loads((HERE / "result.json").read_bytes())["three_jet_traces"]
    total = 0
    for expression, coefficient in traces:
        factors = expression.split()
        slots = [
            (i, j)
            for i, w in enumerate(factors)
            for j, letter in enumerate(w)
            if letter == "A"
        ]
        assert len(slots) == 5
        for removed in combinations(slots, 3):
            total += Fraction(coefficient) * prod(
                tr(["U" if (i, j) in removed else letter for j, letter in enumerate(w)])
                for i, w in enumerate(factors)
            )
    assert (
        total
        == Fraction(32, 5)
        * tr("AAU")
        * (
            tr("UB") ** 2 * tr("CC")
            - 2 * tr("UB") * tr("UC") * tr("BC")
            + tr("UC") ** 2 * tr("BB")
        )
        == Fraction(384, 5)
    )
    return total


def verify(all_colours=False):
    raw = gzip.decompress(SOURCE.read_bytes())
    assert (
        sha256(raw).hexdigest()
        == "ca59e95ecf2ca092722cd612a4cbf51d562f26f06e58a6a7c9d362cdf126ed8a"
    )
    terms = json.loads(raw)["terms"]
    assert len(terms) == 6139
    basis = []
    for i in range(3):
        for j in range(3):
            if i == j and i == 2:
                continue
            matrix = [[0] * 3 for _ in range(3)]
            matrix[i][j] = 1
            if i == j:
                matrix[2][2] = -1
            basis.append(matrix)
    upper = basis[1]
    assert upper[0][1] == 1
    values = response(terms, upper)
    cubic = {m: v for m, v in values.items() if m.bit_count() == 9}
    assert cubic == response(terms, upper, split=True)
    assert cubic == {m: -25477120 for m in (0x7BE, 0xB7E)} | {
        m: -5914624 for m in (0x7EB, 0xBE7, 0xE7B, 0xEB7)
    }
    # Degree eleven selects exactly two removed scalars. Test both
    # extremal adjoint directions for the unique allowed spin-5/2 coupling.
    lower = basis[3]
    assert lower[1][0] == 1
    assert not any(
        m.bit_count() == 11 for v in (values, response(terms, lower)) for m in v
    )
    if all_colours:
        for a, b in combinations_with_replacement(basis, 2):
            direction = [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]
            assert not any(m.bit_count() == 11 for m in response(terms, direction))
        # Homogeneous cubics on eight variables are determined by these
        # 8 + 2*28 + 56 evaluations, in order of the size of their support.
        for a, b, c in combinations_with_replacement(basis, 3):
            direction = [
                [a[i][j] + b[i][j] + c[i][j] for j in range(3)] for i in range(3)
            ]
            result = response(terms, direction)
            assert all(
                result.get(m, 0) == cubic[m] * direction[0][1] ** 3
                for m in (0x7BE, 0x7EB)
            )
    radial = radial_cubic()
    return dict(
        source_terms=len(terms),
        degree_eleven_nonzero=0,
        cubic_ratios={hex(m): str(v) for m, v in sorted(cubic.items())},
        ratio_of_components=str(cubic[0x7BE] / cubic[0x7EB]),
        radial_coefficient=str(radial),
        three_scalar_Wick_colour_ratios={
            hex(m): str(6 * cubic[m] * radial) for m in (0x7BE, 0x7EB)
        },
        colour_pairs_tested=36 if all_colours else 0,
        colour_triples_tested=120 if all_colours else 0,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all-colours", action="store_true")
    if not __debug__:
        parser.error("run with assertions enabled")
    print(json.dumps(verify(parser.parse_args().all_colours), indent=2))
