"""Test a colour-invariant three-fermion, one-scalar contact on the separating Cartan functional."""

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import permutations, product
import json

from cartan_contacts import SCALAR_UNITS, colour_entry, multiply
from cubic_colours import source
from four_fermion_contacts import determinant, projected_basis
from source_fluctuations import add, sign

UNIT = {(0, 0): 1}
INSERTED = 1 << 24
U, V = (1, 0, -1), (0, 1, -1)


@cache
def diagonal_trace(word):
    total = Counter()
    for i in range(3):
        value = UNIT
        for letter in word:
            value = multiply(value, colour_entry(letter, i, i))
        total.update(value)
    return {key: v for key, v in total.items() if v}


@cache
def marked_trace(word, matter, scalar, p, q):
    # The root legs must close in one trace. The fermion uses E_pq and
    # the scalar E_qp. INSERTED is the odd displacement of the fermion.
    terms = []
    for first, psi in enumerate(word):
        if psi != 3 + matter:
            continue
        for second, phi in enumerate(word):
            if phi != scalar:
                continue
            colour = p if first < second else q
            value = UNIT
            for i, letter in enumerate(word):
                if i == first:
                    assert colour == p
                    colour = q
                    value = multiply(value, {(0, INSERTED): 1})
                elif i == second:
                    assert colour == q
                    colour = p
                else:
                    value = multiply(value, colour_entry(letter, colour, colour))
            terms.append(value)
    return add(*terms)


def root_pair(terms, matter, scalar, p, q):
    total = Counter()
    for coefficient, traces in terms:
        words = [tuple(tr) for tr in traces]
        for j, word in enumerate(words):
            marked = marked_trace(word, matter, scalar, p, q)
            if not marked:
                continue
            term = {(0, 0): coefficient}
            for i, tr in enumerate(words):
                term = multiply(term, marked if i == j else diagonal_trace(tr))
                if not term:
                    break
            total.update(term)
    remaining = {}
    for (s, mask), v in total.items():
        if v:
            assert mask & INSERTED and mask.bit_count() == 5
            remaining[s, mask ^ INSERTED] = v
    return remaining


def contact_image(terms):
    result = Counter()
    for p, q in permutations(range(3), 2):
        for i, scalar in product(range(3), repeat=2):
            pair = root_pair(terms, i, scalar, p, q)
            if not pair:
                continue
            remaining = sorted({0, 1, 2} - {i})
            for j, k in permutations(remaining):
                for a, b in product(range(2), repeat=2):
                    weight = (U[q] - U[p]) * V[a] - (V[q] - V[p]) * U[a]
                    if not weight:
                        continue
                    first = 1 << (8 * j + 4 * a)
                    second = 1 << (8 * k + 4 * b)
                    # Move the derivative of flavour i to the front, then
                    # apply the other two in the indicated order.
                    orientation = (-1) ** (i + int(j > k))
                    for (s, mask), v in pair.items():
                        if first & mask and second & mask:
                            tail = mask ^ first
                            parity = sign(first, tail) * sign(second, tail ^ second)
                            result[s + SCALAR_UNITS[2 * scalar + b], tail ^ second] += (
                                weight * orientation * parity * v
                            )
    return {key: v for key, v in result.items() if v}


def independent_specialization(terms):
    from source_fluctuations import matrix_product, trace, wedge

    # The five untouched scalar powers force flavour one. The unsplit
    # matrix substitution below does not use the marked-trace calculation.
    a = (2, -3, 1)
    part = [
        (c, words)
        for c, words in terms
        if not any(letter in (1, 2) for tr in words for letter in tr)
    ]
    assert len(part) == 49
    total = Fraction()
    for p, q in permutations(range(3), 2):
        for i, j, k in permutations(range(3)):
            for weight, direction in ((U[q] - U[p], V), (V[p] - V[q], U)):
                if not weight:
                    continue
                fields = [[[{} for _ in range(3)] for _ in range(3)] for _ in range(6)]
                fields[0][q][p] = {0: 1}
                fields[3 + i][p][q] = {1 << (2 + i): 1}
                for d in range(3):
                    fields[0][d][d] = {0: a[d]}
                    fields[3][d][d] = {1: U[d], 2: V[d]}
                    for h, v in ((j, direction[d]), (k, a[d])):
                        fields[3 + h][d][d][1 << (2 + h)] = v

                @cache
                def word_matrix(w):
                    if len(w) == 1:
                        return fields[w[0]]
                    return matrix_product(word_matrix(w[:-1]), fields[w[-1]])

                @cache
                def traced(w):
                    return trace(word_matrix(w))

                for c, words in part:
                    value = {0: c}
                    for tr in words:
                        value = wedge(value, traced(tuple(tr)))
                    total += weight * value.get(31, 0)
    return total


def verify():
    _, pure, target, functional = projected_basis()
    probes = sorted(functional)
    terms = [(Fraction(c), w) for c, w in source()]
    image = contact_image(terms)
    projected = sum(c * image.get(key, 0) for key, c in functional.items())
    assert len(image) == 120 and projected == Fraction(702, 37)
    sample = sum(
        v * 2 ** (s % 6) * (-3) ** (s // 6 % 6)
        for (s, m), v in image.items()
        if s < 36 and m == 17
    )
    assert sample and independent_specialization(terms) == sample
    columns = pure + [image]
    matrix = [[v.get(key, 0) for v in columns] for key in probes]
    det = determinant(matrix)
    assert det
    weights = []
    for j in range(4):
        replaced = [
            row[:j] + [target.get(key, 0)] + row[j + 1 :]
            for row, key in zip(matrix, probes)
        ]
        weights.append(Fraction(determinant(replaced), det))
    assert weights[-1] * projected == 1
    coordinates = set(target).union(*columns)
    assert all(
        sum(c * v.get(key, 0) for c, v in zip(weights, columns)) == target.get(key, 0)
        for key in coordinates
    )
    return dict(
        cartan_coordinates=len(image),
        distinguished_coefficients=[str(image.get(key, 0)) for key in probes],
        value_on_separator=str(projected),
        independent_five_scalar_component=str(sample),
        mixed_basis_minor=str(det),
        coefficients_reconstructing_radial_component=list(map(str, weights)),
        coordinates_checked=len(coordinates),
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
