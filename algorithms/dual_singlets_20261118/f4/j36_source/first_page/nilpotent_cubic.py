"""Extend the lowest F4 source symbols through the cubic local commuting equations."""

import argparse
from collections import defaultdict
from itertools import combinations_with_replacement
import json
from pathlib import Path

from flint import fmpq, fmpq_mat
from nilpotent_quadratic import (
    add,
    add_term,
    array,
    calculate as quadratic,
    pivots,
    simple_coordinates,
    symbol_monomials,
)
from nilpotent_stabilizer import reflection_word, source_action


def raw_bracket(c, left, right):
    result = {}
    for (axes, i), x in left.items():
        for (other, j), y in right.items():
            if axes[0] in other:
                continue
            key = tuple(sorted((*axes, *other)))
            sign = (-1) ** sum(axes[0] > a for a in other)
            for m, z in c.br[i][j].items():
                add_term(result, (key, m), sign * x * y * z)
    return result


def split(frames, weight, form, closed=False):
    if not form:
        return {}, {}
    basis, inverse, orig, im, number, offset = frames.splittings[weight]
    coefficients = inverse * array(len(basis), [[form.get(e, 0) for e in basis]])
    n = len(im)
    if closed:
        assert not any(coefficients.entries()[n + number :])
    h = {orig[im[i]]: coefficients[i, 0] for i in range(n) if coefficients[i, 0]}
    projection = {
        offset + i: coefficients[n + i, 0]
        for i in range(number)
        if coefficients[n + i, 0]
    }
    return h, projection


def transferred_equations(frames):
    ws = frames.forms[1]
    nm = len(frames.forms[2])
    quad = [{} for _ in range(nm)]
    cubic = [defaultdict(fmpq) for _ in range(nm)]
    lifted = {}
    for i, j in combinations_with_replacement(range(len(ws)), 2):
        w = add(ws[i][0], ws[j][0])
        raw = raw_bracket(frames.colour, ws[i][1], ws[j][1])
        raw = {key: v / (1 + int(i == j)) for key, v in raw.items()}
        h, projection = split(frames, w, raw, True)
        assert frames.bracket(i, 1, j) == {
            m: v * (1 + int(i == j)) for m, v in projection.items()
        }
        if h:
            lifted[i, j] = (w, {key: -v for key, v in h.items()})
        for m, v in projection.items():
            quad[m][i, j] = v
    for i, (wi, form) in enumerate(ws):
        for (j, k), (wjk, h) in lifted.items():
            raw = raw_bracket(frames.colour, form, h)
            _, projection = split(frames, add(wi, wjk), raw)
            monomial = tuple(sorted((i, j, k)))
            for m, v in projection.items():
                add_term(cubic[m], monomial, v)
    # Verify the eliminated two-form equations on three deterministic field configurations.
    # This uses actual nonlinear colour curvatures on the four scalar axes.
    from random import Random

    rng = Random(42197)
    controls = []
    for _ in range(3):
        values = [rng.randint(-2, 2) for _ in ws]
        u = {}
        corr = {}
        for t, (_, form) in zip(values, ws):
            for key, v in form.items():
                add_term(u, key, t * v)
        for (i, j), (_, h) in lifted.items():
            for key, v in h.items():
                add_term(corr, key, values[i] * values[j] * v)
        two = raw_bracket(frames.colour, u, u)
        two = {key: v / 2 for key, v in two.items()}
        three = raw_bracket(frames.colour, u, corr)

        def project_all(raw, closed):
            pieces = defaultdict(dict)
            lookup = {e: w for w, s in frames.splittings.items() for e in s[0]}
            for key, v in raw.items():
                pieces[lookup[key]][key] = v
            result = {}
            preimage = {}
            for w, form in pieces.items():
                h, coeffs = split(frames, w, form, closed)
                result.update(coeffs)
                preimage.update(h)
            return preimage, result

        h, p2 = project_all(two, True)
        _, p3 = project_all(three, False)
        assert h == {key: -v for key, v in corr.items()}
        for m in range(nm):
            assert p2.get(m, 0) == sum(
                v * values[i] * values[j] for (i, j), v in quad[m].items()
            )
            assert p3.get(m, 0) == sum(
                v * values[i] * values[j] * values[k]
                for (i, j, k), v in cubic[m].items()
            )
        controls.append(dict(quadratic=len(p2), cubic=len(p3)))
    return quad, cubic, controls


def multiply_symbols(symbols, equations):
    return [
        {tuple(sorted((*f, *term))): v for term, v in equations[m].items()}
        for f, m in symbols
    ]


def kernel(matrix):
    red, columns = pivots(matrix)
    vectors = []
    for j in sorted(set(range(matrix.ncols())) - set(columns)):
        v = [fmpq(0)] * matrix.ncols()
        v[j] = 1
        for i, k in enumerate(columns):
            v[k] = -red[i, j]
        vectors.append(v)
    return array(matrix.ncols(), vectors)


def calculate(q, short, order):
    base, frames, mons, inc, out = quadratic(q, short, order, True)
    roots = simple_coordinates(frames.colour)
    target = tuple(
        4 * s - 2 * r for s, r in zip(add(*(roots[i] for i in q)), roots[q[short]])
    )
    next_mons = symbol_monomials(frames, target, order + 1)
    quad, cubic, controls = transferred_equations(frames)
    scores = multiply_symbols(mons, cubic)
    corrections = multiply_symbols(next_mons, quad)
    coords = sorted({u for col in scores + corrections for u in col})
    build = lambda columns: array(
        len(coords), [[col.get(k, 0) for k in coords] for col in columns]
    )
    score, correction = build(scores), build(corrections)
    cycles = kernel(out)

    def span(m, n):
        assert m.nrows() == n.nrows()
        return array(
            m.nrows(),
            [
                [a[i, j] for i in range(a.nrows())]
                for a in (m, n)
                for j in range(a.ncols())
            ],
        )

    corr_rank = correction.rank()
    assert span(correction, score * inc).rank() == corr_rank
    # The involution preserves the quadratic differential. Finite averaging
    # supplies invariant representatives for all invariant quotient classes.
    averaged = cycles
    for word in reflection_word(frames.colour, roots, q, short):
        if word:
            a, char, receipt = source_action(q, short, frames, mons, word)
            identity = fmpq_mat(
                [[int(i == j) for j in range(len(mons))] for i in range(len(mons))]
            )
            averaged = averaged * kernel((a - char * identity) * averaged)
    n = span(inc, averaged).rank() - inc.rank()
    rank = span(correction, score * averaged).rank() - corr_rank
    post = span(correction, score * cycles).rank() - corr_rank
    return dict(
        root_colours=q,
        short_axis=short,
        order=order,
        initial_symbols=len(mons),
        first_corrections=len(next_mons),
        all_marked_quadratic_monomials=sum(map(len, quad)),
        all_marked_cubic_monomials=sum(map(len, cubic)),
        direct_curvature_controls=controls,
        target_monomials=len(coords),
        quadratic_boundary_rank=corr_rank,
        cubic_score_nonzero_columns=sum(map(bool, scores)),
        ordinary_obstruction_rank=post,
        finite_stabilizer_before=n,
        finite_stabilizer_obstruction_rank=rank,
        finite_stabilizer_survivor_bound=n - rank,
        correction_matrix=[[str(c.get(k, 0)) for c in corrections] for k in coords],
        cubic_score_matrix=[[str(c.get(k, 0)) for c in scores] for k in coords],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    data = calculate((8, 19, 27, 47), 1, 6)
    a.out.write_text(json.dumps(data, indent=2) + "\n")
    print(
        data["finite_stabilizer_before"],
        data["finite_stabilizer_obstruction_rank"],
        data["finite_stabilizer_survivor_bound"],
    )
