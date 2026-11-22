"""Check the scalar-discriminant homotopy on the pinned SU(3) R-007 operator."""

from collections import Counter
from fractions import Fraction
from functools import cache, lru_cache
from hashlib import sha256
from math import comb, prod
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import json
import runpy

from cartan_contacts import cartan_scalar_expansion
from cubic_colours import source
from source_fluctuations import add, scale
from verify import A, B, C, HASHES, blob, endpoint


def odd(slot):
    # Three components (diagonal, upper, lower) for each of B,C,psi1,psi2,psi3,f.
    return 3 <= slot // 3 <= 5


@lru_cache(65536)
def order(word):
    exterior = [x for x in word if odd(x)]
    if len(exterior) != len(set(exterior)):
        return (), 0
    swaps = sum(x > y for i, x in enumerate(exterior) for y in exterior[i + 1 :])
    return tuple(sorted(word)), (-1) ** swaps


def multiply(left, right):
    result = Counter()
    for w, a in left.items():
        for v, b in right.items():
            key, sign = order(w + v)
            if sign:
                result[key] += a * b * sign
    return {k: v for k, v in result.items() if v}


def act(terms, rules, divide_roots=False):
    result = Counter()
    for word, value in terms.items():
        if divide_roots:
            count = sum(bool(x % 3) for x in word)
            if not count:
                continue
            value = Fraction(value, count)
        for i, letter in enumerate(word):
            for replacement, weight in rules.get(letter, {}).items():
                key, sign = order(word[:i] + replacement + word[i + 1 :])
                if sign:
                    result[key] += value * sign * weight
            if odd(letter):
                value = -value
    return {k: v for k, v in result.items() if v}


def entries(letter):
    diagonal, upper, lower = [(3 * letter + i,) for i in range(3)]
    return [[{diagonal: 1}, {upper: 1}], [{lower: 1}, {diagonal: -1}]]


def matrix_product(left, right):
    return [
        [add(*(multiply(left[i][k], right[k][j]) for k in range(2))) for j in range(2)]
        for i in range(2)
    ]


def restricted_source(terms, eigenvalues):
    matrices = [[[{(): eigenvalues[0]}, {}], [{}, {(): eigenvalues[1]}]]]
    matrices.extend(entries(i) for i in range(1, 6))

    @cache
    def matrix(word):
        if len(word) == 1:
            return matrices[word[0]]
        return matrix_product(matrix(word[:-1]), matrices[word[-1]])

    @cache
    def trace(word):
        m = matrix(word)
        result = add(m[0][0], m[1][1])
        if all(i == 0 for i in word):
            result = add(result, {(): eigenvalues[2] ** len(word)})
        return result

    result = Counter()
    for coefficient, traces in terms:
        polynomial = {(): Fraction(coefficient)}
        for word in traces:
            polynomial = multiply(polynomial, trace(tuple(word)))
        result.update(polynomial)
    return {w: v for w, v in result.items() if v}


def split_charge(kernel, eigenvalues):
    # Restore the middle public fermion sign before reading the commutators.
    masks = [(0,), (1,), (2,), (1, 2), (0, 2), (0, 1), (0, 1, 2)]
    matrices = [[[{(): eigenvalues[0]}, {}], [{}, {(): eigenvalues[1]}]]]
    matrices.extend(entries(i) for i in range(1, 7))
    linear, quadratic = {}, {}
    for i in range(1, 7):
        components = [[Counter(), Counter()], [Counter(), Counter()]]
        for c, pair in kernel((0, 0, masks[i])):
            left, right = (masks.index(x[2]) for x in pair)
            coeff = c * (-1) ** ((i == 4) + (left == 4) + (right == 4))
            result = matrix_product(matrices[left], matrices[right])
            for row in range(2):
                for col in range(2):
                    components[row][col].update(scale(result[row][col], coeff))
        for axis, (row, col) in enumerate(((0, 0), (0, 1), (1, 0))):
            component = components[row][col]
            slot = 3 * i + axis
            linear[slot] = {w: c for w, c in component.items() if c and len(w) == 1}
            quadratic[slot] = {w: c for w, c in component.items() if c and len(w) == 2}
            assert all(len(w) in (1, 2) for w, c in component.items() if c)
        assert not add(components[0][0], components[1][1])
    h = {}
    for preimage, image in linear.items():
        if image:
            (((letter,), coefficient),) = image.items()
            h[letter] = {(preimage,): Fraction(1, coefficient)}
    assert len(h) == 6
    for slot in linear:
        x = {(slot,): 1}
        anticommutator = add(act(act(x, linear), h), act(act(x, h), linear))
        assert anticommutator == (x if slot % 3 else {})
    return linear, quadratic, h


def check_primitive(terms, kernel, eigenvalues):
    linear, quadratic, h = split_charge(kernel, eigenvalues)
    full = {i: add(linear[i], quadratic[i]) for i in linear}
    assert all(not act(v, full) for v in full.values())
    tree = restricted_source(terms, eigenvalues)
    assert tree and not act(tree, full)
    residual, primitive = tree.copy(), {}
    steps = []
    for degree in range(5, 9):
        row = {w: c for w, c in residual.items() if len(w) == degree}
        assert not act(row, linear) and all(any(x % 3 for x in w) for w in row)
        part = act(row, h, divide_roots=True)
        assert act(part, linear) == row
        primitive.update(part)
        residual = add(residual, scale(act(part, full), -1))
        steps.append([degree, len(row), len(part)])
    assert not residual and act(primitive, full) == tree
    return dict(
        eigenvalues=eigenvalues,
        nonzero_source_coefficients=len(tree),
        nonzero_primitive_coefficients=len(primitive),
        degree_residual_primitive=steps,
        final_residual=0,
    )


def finite_difference(order_bound, exponent):
    # Lagrange evaluation at zero from t=exponent,...,exponent+order_bound.
    weights = [
        (-1) ** i
        * comb(exponent + i - 1, i)
        * comb(exponent + order_bound, order_bound - i)
        for i in range(order_bound + 1)
    ]
    for j in range(order_bound + 1):
        assert sum(c * comb(exponent + i, j) for i, c in enumerate(weights)) == int(
            j == 0
        )
    return weights


def verify():
    terms = source()
    counts = Counter(sum(x != 0 for tr in words for x in tr) for _, words in terms)
    assert counts == {5: 49, 6: 866, 7: 3096, 8: 2128}
    assert all(
        sum((0, 1, 1, 2, 1, 1)[x] for tr in words for x in tr) == 8
        for _, words in terms
    )
    _, root_support = cartan_scalar_expansion(terms)
    assert root_support == [0, 0, 3986, 0, 0, 0]
    data = {name: blob("dot-r007-t2-primary/" + name) for name in HASHES}
    assert all(sha256(v).hexdigest() == HASHES[k] for k, v in data.items())
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "replay.py"
        path.write_bytes(data["replay.py"])
        engine = runpy.run_path(str(path))
        outgoing = endpoint(SimpleNamespace(**engine), data)
    values = {A: (1, 2, -3), B: (1, 0, -1), C: (0, 1, -1)}
    selected = {
        w: c for w, c in outgoing.items() if all(x in values for tr in w for x in tr)
    }
    value = sum(
        c * prod(sum(prod(values[x][i] for x in tr) for i in range(3)) for tr in w)
        for w, c in selected.items()
    )
    assert len(outgoing) == 5080 and len(selected) == 5 and value == Fraction(-14112, 5)
    roots = values[A]
    delta = prod((x - y) ** 2 for i, x in enumerate(roots) for y in roots[i + 1 :])
    assert delta == 400
    estimates = [
        [n, n - min(counts) + 1, n + (n - min(counts) + 1)] for n in range(5, 9)
    ]
    bound = (max(v[-1] for v in estimates) + 1) // 2
    assert bound == 6
    return dict(
        source_words_by_nonbackground_degree=dict(sorted(counts.items())),
        diagonal_scalar_terms_by_off_diagonal_fermion_count=root_support,
        degree_inverse_gaps_projector_plus_gap_bound=estimates,
        sufficient_discriminant_power=bound,
        regular_cartan_discriminant=delta,
        regular_target_value=str(value),
        ward_weights_at_orders_zero_through_four=[
            finite_difference(m, bound) for m in range(5)
        ],
        block_checks=[
            check_primitive(terms, engine["q_letter"], v) for v in (roots, (1, -3, 2))
        ],
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
