"""Finite Weyl stabilizer on the lowest transverse F4 source symbols."""

from collections import defaultdict, deque
from functools import cache
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path
import argparse

from flint import fmpq, fmpq_mat
from nilpotent_flag_orders import reflections
from nilpotent_quadratic import (
    add_term,
    array,
    calculate as quadratic,
    pivots,
    simple_coordinates,
)


def acted(matrix, vector):
    res = {}
    for source, x in vector.items():
        for target, y in matrix[source].items():
            add_term(res, target, x * y)
    return res


def reflection_word(colour, roots, quartet, short):
    basis = tuple(range(4, 52))
    action = reflections(colour, roots)
    words = {basis: ()}
    todo = deque([basis])
    while todo:
        w = todo.popleft()
        for i, s in enumerate(action):
            out = tuple(s[j] for j in w)
            if out not in words:
                words[out] = (*words[w], i)
                todo.append(out)
    assert len(words) == 1152
    return [
        word
        for w, word in words.items()
        if w[quartet[short] - 4] == quartet[short]
        and {w[i - 4] for i in quartet} == set(quartet)
    ]


def native_lift(c, word):
    @cache
    def simple(i):
        e, f = c.data["raising"][i], c.data["lowering"][i]
        h = c.br[e][f]
        coefficient = sum(v * c.weights[e][j] for j, v in h.items())

        def exp_action(root, scale, initial):
            result = initial.copy()
            term = initial
            degree = 1
            while term:
                term = {
                    k: (scale * x) / degree for k, x in acted(c.br[root], term).items()
                }
                for k, v in term.items():
                    add_term(result, k, v)
                degree += 1
                assert degree < 10
            return result

        def apply(v):
            v = exp_action(e, fmpq(1), v)
            v = exp_action(f, -fmpq(2, coefficient), v)
            return exp_action(e, fmpq(1), v)

        return [apply({j: fmpq(1)}) for j in range(52)]

    identity = [{j: fmpq(1)} for j in range(52)]
    lift = identity
    for i in word:
        lift = [acted(simple(i), v) for v in lift]
    for i, j in combinations_with_replacement(range(52), 2):
        check = {}
        for a, x in lift[i].items():
            for b, y in lift[j].items():
                for k, val in c.br[a][b].items():
                    add_term(check, k, x * y * val)
        assert check == acted(lift, c.br[i][j])
    power = identity
    for exponent in range(1, 97):
        power = [acted(lift, v) for v in power]
        if power == identity:
            return lift, exponent
    raise ValueError("native lift failed the finite order test")


def source_action(q, short, frames, monomials, word):
    c = frames.colour
    roots = simple_coordinates(c)
    lift, exponent = native_lift(c, word)
    assert all(len(lift[i]) == 1 for i in q)
    images = [next(iter(lift[i].items())) for i in q]
    perm = [q.index(j) for j, _ in images]
    scales = [x for _, x in images]
    # The two full columns and two columns on the complementary three-plane
    # give the fourth power of a long-axis permutation sign.
    character = 1 / prod(x ** (2 if i == short else 4) for i, x in enumerate(scales))
    assert character**exponent == 1
    reverse = {}
    for degree in (1, 2):
        images = defaultdict(dict)
        for i, (_, form) in enumerate(frames.forms[degree]):
            t = {}
            for (axes, colour), v in form.items():
                mapped = [perm[a] for a in axes]
                sign = (-1) ** sum(
                    a > b for j, a in enumerate(mapped) for b in mapped[j + 1 :]
                )
                value = sign * v / prod(scales[a] for a in axes)
                for k, x in lift[colour].items():
                    add_term(t, (tuple(sorted(mapped)), k), value * x)
            weights = {
                tuple(
                    r - sum(roots[q[a]][j] for a in axes)
                    for j, r in enumerate(roots[k])
                )
                for axes, k in t
            }
            assert len(weights) == 1
            for j, x in frames.project(degree, weights.pop(), t).items():
                images[j][i] = x
        reverse[degree] = images
    coord = {m: i for i, m in enumerate(monomials)}
    matrix = fmpq_mat(len(coord), len(coord))
    for j, (factors, mark) in enumerate(monomials):
        values = {(): fmpq(1)}
        for i in factors:
            out = {}
            for key, v in values.items():
                for t, x in reverse[1][i].items():
                    add_term(out, tuple(sorted((*key, t))), v * x)
            values = out
        for a, x in values.items():
            for m, y in reverse[2][mark].items():
                matrix[coord[a, m], j] += x * y
    identity = fmpq_mat(
        [[int(i == j) for j in range(len(coord))] for i in range(len(coord))]
    )
    assert matrix**exponent == identity
    return (
        matrix,
        character,
        dict(
            word=word,
            order_on_colour=exponent,
            axis_permutation=perm,
            axis_scalings=list(map(str, scales)),
            character=str(character),
        ),
    )


def calculate(q, short, order):
    base, frames, mons, inc, out = quadratic(q, short, order, True)
    words = reflection_word(frames.colour, simple_coordinates(frames.colour), q, short)
    n = len(mons)
    _, independent = pivots(inc)
    images = [[inc[i, j] for i in range(n)] for j in independent]
    red, cols = pivots(out)
    ker = []
    for j in sorted(set(range(n)) - set(cols)):
        v = [fmpq(0)] * n
        v[j] = 1
        for i, k in enumerate(cols):
            v[k] = -red[i, j]
        ker.append(v)
    all_cols = images + ker
    _, keep = pivots(array(n, all_cols))
    basis = array(n, [all_cols[j] for j in keep])
    _, selected = pivots(basis.transpose())
    square = fmpq_mat([[basis[i, j] for j in range(len(keep))] for i in selected])
    quotient = base["surviving_dimension"]
    constraints = []
    records = []
    for w in words:
        if not w:
            continue
        a, char, receipt = source_action(q, short, frames, mons, w)
        transformed = a * basis
        induced = square.inv() * fmpq_mat(
            [[transformed[i, j] for j in range(len(keep))] for i in selected]
        )
        assert basis * induced == transformed
        m = len(images)
        assert all(not induced[i, j] for i in range(m, len(keep)) for j in range(m))
        block = [
            [induced[m + i, m + j] for j in range(quotient)] for i in range(quotient)
        ]
        for i, row in enumerate(block):
            constraints.append([v - char * int(i == j) for j, v in enumerate(row)])
        receipt["action_on_quadratic_cohomology"] = [
            [str(v) for v in row] for row in block
        ]
        records.append(receipt)
    rank = fmpq_mat(constraints).rank() if constraints else 0
    return dict(
        quadratic=base,
        stabilizing_weyl_elements=len(words),
        nonidentity_actions=records,
        finite_stabilizer_dimension=quotient - rank,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quartet", nargs=4, type=int, default=(8, 19, 27, 47))
    parser.add_argument("--short-axis", type=int, choices=range(4), default=1)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--out", type=Path, required=True)
    a = parser.parse_args()
    if not __debug__ or not 0 <= a.order <= 8:
        parser.error("enable assertions and specify an order from zero through eight")
    data = calculate(tuple(a.quartet), a.short_axis, a.order)
    a.out.write_text(json.dumps(data, indent=2) + "\n")
    print(data["quadratic"]["surviving_dimension"], data["finite_stabilizer_dimension"])
