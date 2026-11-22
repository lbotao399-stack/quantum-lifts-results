"""Compute the top first pairing and the six-substitution orbit in SO(19)."""

# ruff: noqa: E402 -- load the existing trace runtime from its repository location.
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, product
import json
from math import comb, factorial, prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "scalar_multipliers"))
from verify import GHOSTS, PHI, PfaffianCoefficient, encode, pullback, so_canon, words
from first_family import coefficient
from raising import raw_word, substitute, terminal
from q1d_exact import pair_terms
from so8_level24_descendant import check_letter, substitution
from orthogonal_glueball_relation_span import determinant
from verify import unit

STEPS = tuple(product(range(2), range(3)))


class PairedPaths:
    """Specialize the trace logarithm to unlabelled copies of a single double path."""

    def __init__(self, power):
        self.doubles = power - 4
        self.projector = PfaffianCoefficient(3, self.doubles)

    @cache
    def trace_coefficient(self, word):
        result = Counter()
        orbit = {
            seq[i:] + seq[:i] for seq in (word, word[::-1]) for i in range(len(seq))
        }
        for seq in orbit:
            orientation, canon = so_canon((seq,))
            assert orientation and canon == (word,)
            letters = tuple((PHI + GHOSTS).index(t) for t in seq)
            axes = []
            doubles = 0
            i = 0
            while i < len(letters):
                if letters[i] < 3:
                    axes.append(letters[i])
                    i += 1
                    continue
                pair = letters[i : i + 2]
                if pair not in ((3, 4), (4, 3)):
                    break
                orientation *= 1 if pair == (4, 3) else -1
                doubles += 1
                i += 2
            else:
                blocks = len(axes) + doubles
                multiplier = Fraction(orientation * (-1) ** (blocks + 1), 2 * blocks)
                for (mask, degree), value in self.projector.aux_trace(
                    tuple(axes)
                ).items():
                    assert mask == 0
                    result[degree, doubles] += multiplier * value
        return {key: value for key, value in result.items() if value}

    @cache
    def product_coefficient(self, traces):
        if not traces:
            return {(0, 0): 1}
        result = Counter()
        for (i, j), c in self.product_coefficient(traces[:-1]).items():
            for (a, b), v in self.trace_coefficient(traces[-1]).items():
                if i + a in self.projector.divisors and j + b <= self.doubles:
                    result[i + a, j + b] += c * v
        return {key: value for key, value in result.items() if value}

    def coefficient(self, traces):
        multiplicities = prod(factorial(n) for n in Counter(traces).values())
        value = sum(
            c * self.projector.weights.get(i, 0)
            for (i, j), c in self.product_coefficient(traces).items()
            if j == self.doubles
        )
        return Fraction(factorial(self.doubles), multiplicities) * value


def check_actions():
    jets = {
        n: [
            (s, t, f)
            for s in range(n + 1)
            for t in range(n - s + 1)
            for f in combinations(range(3), n - s - t)
        ]
        for n in (1, 2)
    }
    letters = jets[1] + jets[2]
    pairs = (
        list(product(jets[1], jets[1]))
        + list(product(jets[1], jets[2]))
        + list(product(jets[2], jets[1]))
    )
    for step in STEPS:
        for x in letters:
            check_letter(x, *step, True)
            for other in STEPS:
                image = Counter()
                for first, second in ((step, other), (other, step)):
                    c, y = substitution(x, *first, True)
                    if c:
                        v, z = substitution(y, *second, True)
                        image[z] += c * v
                assert not any(image.values())
        for x, y in pairs:
            for pos, letter in enumerate((x, y)):
                c, z = substitution(letter, *step, True)
                if c:
                    assert not (pair_terms(z, y) if pos == 0 else pair_terms(x, z))
            for _, left, right in pair_terms(x, y):
                assert not substitution(left, *step, True)[0]
                assert not substitution(right, *step, True)[0]
    return len(letters), len(pairs)


def top_pairing(colours, power, direct=False):
    traces = {}
    for monomial in words((power - 3, power - 3, 2, 2, 2)):
        if any(sum(factor) != 2 for factor in monomial):
            continue
        original = raw_word(monomial)
        ((word, orientation),) = original.items()
        for step in STEPS:
            original = substitute(original, *step, True)
        value = sum(c * terminal(w)[0] for w, c in original.items())
        if value:
            traces[word] = orientation * value
    dual, report = pullback(colours, traces, 1)
    assert power + 5 > colours and power < colours - 2
    paths = PairedPaths(power)
    if direct:
        n = power - 4
        spec = dict(
            paths=["mp"] * n,
            matrices=[unit(4, 3)] * n,
            scalar_directions=[[]] * n,
            single=((), ()),
        )
        native = determinant(power, spec)
        native.project, native.letters, native.pf_degree = (
            PfaffianCoefficient(3, 0),
            2 * power + 1,
            power + 5,
        )
        for w in dual:
            assert paths.coefficient(w) == native.coefficient(encode(w))
        native.clear()
    value = sum(v * paths.coefficient(w) for w, v in dual.items())
    assert value == (-1) ** power * 1440 * factorial(power - 4) * coefficient(
        colours, power
    )
    paths.trace_coefficient.cache_clear()
    paths.product_coefficient.cache_clear()
    return dict(
        colours=colours,
        power=power,
        top=str(value),
        checked_columns=report["checked_columns"],
        ordered_terms=len(dual),
        direct_word_comparison=direct,
    )


def multiplets():
    # Exterior Cauchy. Two-row diagrams inside the 2 by 3 rectangle.
    # The SL(3) Dynkin labels are those for the dual of the transposed diagram.
    blocks = []
    for degree in range(7):
        rows = []
        for a in range(4):
            b = degree - a
            if not 0 <= b <= a:
                continue
            t = [int(a >= j) + int(b >= j) for j in (1, 2, 3)]
            p, q = t[1] - t[2], t[0] - t[1]
            dimension = (a - b + 1) * (p + 1) * (q + 1) * (p + q + 2) // 2
            rows.append(
                dict(
                    twice_dotted_spin=a - b, flavour_Dynkin=[p, q], dimension=dimension
                )
            )
        assert sum(row["dimension"] for row in rows) == comb(6, degree)
        blocks.append(rows)
    assert sum(map(len, blocks)) == 10
    return blocks


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Use assertions")
    check = check_actions()
    cases = [
        top_pairing(n, k, direct=k <= 6)
        for n, k in ((7, 4), (9, 5), (10, 6), (19, 15), (19, 16))
    ]
    print(
        json.dumps(
            dict(
                letters=check[0],
                ordered_pairs=check[1],
                cases=cases,
                highest_weights_by_degree=multiplets(),
            ),
            indent=2,
        )
    )
