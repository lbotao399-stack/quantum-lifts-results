"""Expand the SU(3) bare source at diagonal scalars and test the support of a local contact."""

from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import product
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import json
import runpy

from cubic_colours import source
from source_fluctuations import add, sign, trace
from verify import HASHES, blob, endpoint


def radial_support(terms):
    flavours = set(range(3))
    bare = [(0, 0, (i,)) for i in range(3)]
    bare.extend((0, 0, tuple(sorted(flavours - {i}))) for i in range(3))
    source_contents = {
        tuple(sum(tr.count(i) for tr in words) for i in range(6)) for _, words in terms
    }
    data = {name: blob("dot-r007-t2-primary/" + name) for name in HASHES}
    assert all(
        sha256(value).hexdigest() == HASHES[name] for name, value in data.items()
    )
    with TemporaryDirectory() as temp:
        path = Path(temp) / "replay.py"
        path.write_bytes(data["replay.py"])
        engine = SimpleNamespace(**runpy.run_path(str(path)))
        target = endpoint(engine, data)

    changed_counts = Counter()
    jet_pairs = Counter()
    selected = {}
    for words, coefficient in target.items():
        letters = sum(words, ())
        counts = [letters.count(letter) for letter in bare]
        unchanged = max(
            sum(min(a, b) for a, b in zip(s, counts)) for s in source_contents
        )
        changed = 10 - unchanged
        changed_counts[changed] += 1
        if changed == 3:
            jets = tuple(sorted((i, j) for i, j, theta in letters if not theta))
            assert len(jets) == 2 and all(i + j == 2 for i, j in jets)
            jet_pairs[jets] += 1
            selected[words] = coefficient
    assert sum(changed_counts.values()) == 5080
    diagonal = engine.rho(selected)
    assert changed_counts == {3: 69, 4: 2007, 5: 2877, 6: 126, 7: 1}
    assert len(diagonal) == 240
    return dict(
        bare_contents=len(source_contents),
        target_words_by_minimum_changed_inputs=sorted(changed_counts.items()),
        dotted_jet_pairs=[[jets, n] for jets, n in sorted(jet_pairs.items())],
        selected_target_cartan_monomials=len(diagonal),
    ), diagonal


# Six ordinary scalar coordinates. Base six is lossless since every term has
# scalar degree five. Exterior bits enumerate the eight independent entries
# of each matter fermion in row-major order, omitting the last diagonal.
SCALAR_UNITS = tuple(6**i for i in range(6))
OFF_DIAGONAL = sum(
    1 << (8 * flavour + 3 * row + col)
    for flavour, row, col in product(range(3), repeat=3)
    if row != col
)


def multiply(left, right):
    total = Counter()
    for (s, mask), value in left.items():
        for (t, other), weight in right.items():
            if not mask & other:
                total[s + t, mask | other] += value * weight * sign(mask, other)
    return {key: value for key, value in total.items() if value}


def colour_entry(letter, row, col):
    if letter < 3:
        if row != col:
            return {}
        if row == 2:
            return {(SCALAR_UNITS[2 * letter + j], 0): -1 for j in range(2)}
        return {(SCALAR_UNITS[2 * letter + row], 0): 1}
    start = 8 * (letter - 3)
    if row == col == 2:
        return {(0, 1 << (start + j)): -1 for j in (0, 4)}
    return {(0, 1 << (start + 3 * row + col)): 1}


@cache
def word_matrix(word):
    if len(word) == 1:
        return [[colour_entry(word[0], i, j) for j in range(3)] for i in range(3)]
    left, right = word_matrix(word[:-1]), word_matrix(word[-1:])
    return [
        [add(*(multiply(left[i][k], right[k][j]) for k in range(3))) for j in range(3)]
        for i in range(3)
    ]


@cache
def trace_polynomial(word):
    return trace(word_matrix(word))


def cartan_scalar_expansion(terms):
    denominators = {Fraction(c).denominator for c, _ in terms}
    scale = lcm(*denominators)
    coefficients = Counter()
    for c, words in terms:
        factor = Fraction(c) * scale
        assert factor.denominator == 1
        term = {(0, 0): factor.numerator}
        for tr in words:
            term = multiply(term, trace_polynomial(tuple(tr)))
            if not term:
                break
        coefficients.update(term)
    polynomial = {key: Fraction(v, scale) for key, v in coefficients.items() if v}
    support = Counter((mask & OFF_DIAGONAL).bit_count() for _, mask in polynomial)
    assert all(mask.bit_count() == 5 for _, mask in polynomial)
    assert all(sum(s // unit % 6 for unit in SCALAR_UNITS) == 5 for s, _ in polynomial)
    return polynomial, [support[k] for k in range(6)]


def numeric_check(polynomial):
    # Independently evaluate six scalar coordinates and five common exterior
    # generators by multiplying ordinary Grassmann-valued colour matrices.
    from source_fluctuations import matrix_product, wedge
    from random import Random

    rng = Random(19)
    values = []
    for letter in range(6):
        matrix = [
            [
                [rng.randint(-2, 2) for _ in range(1 if letter < 3 else 5)]
                for _ in range(3)
            ]
            for _ in range(3)
        ]
        matrix[2][2] = [-a - b for a, b in zip(matrix[0][0], matrix[1][1])]
        values.append(matrix)
    scalar = [values[i][j][j][0] for i in range(3) for j in range(2)]
    fermion = [
        {1 << k: x for k, x in enumerate(values[i][j][h]) if x}
        for i in range(3, 6)
        for j, h in product(range(3), repeat=2)
        if (j, h) != (2, 2)
    ]

    @cache
    def exterior(mask):
        if not mask:
            return {0: 1}
        bit = mask & -mask
        return wedge(fermion[bit.bit_length() - 1], exterior(mask ^ bit))

    projected = Counter()
    for (s, mask), value in polynomial.items():
        factor = exterior(mask).get(31, 0) * value
        for v, unit in zip(scalar, SCALAR_UNITS):
            factor *= v ** (s // unit % 6)
        projected[(mask & OFF_DIAGONAL).bit_count()] += factor

    terms = source()
    direct = []
    for parameter in (0, 1, -2):
        fields = []
        for i, matrix in enumerate(values):
            if i < 3:
                field = [
                    [{0: e[0]} if j == h else {} for h, e in enumerate(row)]
                    for j, row in enumerate(matrix)
                ]
            else:
                field = [
                    [
                        {
                            1 << k: v * (1 if j == h else parameter)
                            for k, v in enumerate(e)
                            if v
                        }
                        for h, e in enumerate(row)
                    ]
                    for j, row in enumerate(matrix)
                ]
            fields.append(field)

        @cache
        def ordered(w):
            if len(w) == 1:
                return fields[w[0]]
            return matrix_product(ordered(w[:-1]), fields[w[-1]])

        @cache
        def traced(w):
            return trace(ordered(w))

        total = Fraction()
        for c, words in terms:
            term = {0: Fraction(c)}
            for tr in words:
                term = wedge(term, traced(tuple(tr)))
                if not term:
                    break
            total += term.get(31, 0)
        assert total == sum(v * parameter**k for k, v in projected.items())
        direct.append([parameter, str(total)])
    return direct


def cubic_colour():
    # Coefficients of tr(X^2 U X V - X^2 V X U) for U=H0 and V=H1.
    u, v = (1, 0, -1), (0, 1, -1)
    coefficients = Counter()
    for i, j, k in product(range(3), repeat=3):
        word = tuple(sorted((3 * i + j, 3 * j + k, 3 * k + i)))
        coefficients[word] += u[k] * v[i] - v[k] * u[i]
    result = {w: c for w, c in coefficients.items() if c}
    assert result == {(1, 5, 6): 3, (2, 3, 7): -3}
    return [[list(w), c] for w, c in sorted(result.items())]


def verify():
    terms = source()
    assert len(terms) == 6139
    polynomial, support = cartan_scalar_expansion(terms)
    assert support == [0, 0, 3986, 0, 0, 0]
    radial, _ = radial_support(terms)
    return dict(
        source_words=len(terms),
        source_monomials_by_off_diagonal_fermion_degree=support,
        independent_colour_evaluation=numeric_check(polynomial),
        cubic_colour_cycles=cubic_colour(),
        **radial,
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
