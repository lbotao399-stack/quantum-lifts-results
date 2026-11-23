"""Evaluate local source symbols on exact straight arcs of commuting F4 quadruples."""

import argparse
from collections import defaultdict
from functools import cache
from itertools import combinations, combinations_with_replacement
import json
from math import prod
from pathlib import Path
from random import Random

from flint import fmpq, fmpq_mat
from nilpotent_cubic import kernel, raw_bracket
from nilpotent_quadratic import add_term, array, calculate as quadratic


def linear_combination(forms, values):
    out = {}
    for i, coefficient in values.items():
        for key, x in forms[i][1].items():
            add_term(out, key, coefficient * x)
    return out


def direct_curve(c, q, displacement, mark):
    """Check the actual equations at three points, sufficient for a quadratic polynomial."""

    def comm(x, y):
        v = {}
        for i, a in x.items():
            for j, b in y.items():
                for k, z in c.br[i][j].items():
                    add_term(v, k, a * b * z)
        return v

    for t in (0, 1, -1):
        x = [{j: 1} for j in q]
        for ((axis,), k), v in displacement.items():
            add_term(x[axis], k, t * v)
        assert all(not comm(x[i], x[j]) for i, j in combinations(range(4), 2))
        y = defaultdict(dict)
        for (axes, k), v in mark.items():
            y[axes][k] = v
        for axes in combinations(range(4), 3):
            boundary = {}
            for s, i in enumerate(axes):
                for k, v in comm(x[i], y[axes[:s] + axes[s + 1 :]]).items():
                    add_term(boundary, k, (-1) ** s * v)
            assert not boundary


def calculate(q, short, order, samples=4):
    base, frames, monomials, incoming, outgoing = quadratic(q, short, order, True)
    n = len(monomials)
    marks = sorted({m for f, m in monomials})
    tangent, two = frames.forms[1:3]
    threshold = n - incoming.rank()
    independent = []

    def append(row):
        row = dict(row)
        for pivot, known in independent:
            mult = row.get(pivot, 0)
            if mult:
                for k, v in known.items():
                    add_term(row, k, -mult * v)
        if row:
            pivot = min(row)
            s = row[pivot]
            independent.append((pivot, {k: fmpq(v) / s for k, v in row.items()}))
            return True
        return False

    for i in range(outgoing.nrows()):
        append({j: outgoing[i, j] for j in range(n) if outgoing[i, j]})

    @cache
    def pair(i, j):
        v = raw_bracket(frames.colour, tangent[i][1], tangent[j][1])
        return {k: x / (1 + int(i == j)) for k, x in v.items()}

    @cache
    def marked(i, m):
        return raw_bracket(frames.colour, tangent[i][1], two[m][1])

    by_mark = defaultdict(set)
    for f, m in monomials:
        by_mark[m].add(tuple(sorted(set(f))))
    minimal = set().union(*by_mark.values()) if by_mark else set()
    extra = {
        tuple(sorted(set(f) | set(g)))
        for choices in by_mark.values()
        for f, g in combinations(choices, 2)
    }
    supports = sorted(minimal, key=lambda u: (len(u), u)) + sorted(
        extra - minimal, key=lambda u: (len(u), u)
    )
    recorded = []
    tried = 0
    rng = Random(144362)
    for axes in supports:
        for trial in range(samples):
            if len(independent) == threshold:
                break
            u = {i: (1 if trial == 0 else rng.choice((-2, -1, 1, 2))) for i in axes}
            curvature = {}
            for i, j in combinations_with_replacement(axes, 2):
                for key, v in pair(i, j).items():
                    add_term(curvature, key, u[i] * u[j] * v)
            if curvature:
                continue
            tried += 1
            equations = []
            for m in marks:
                col = {}
                for i, a in u.items():
                    for k, v in marked(i, m).items():
                        add_term(col, k, a * v)
                equations.append(col)
            keys = sorted({k for col in equations for k in col})
            allowed = kernel(
                array(len(keys), [[col.get(k, 0) for k in keys] for col in equations])
            )
            for j in range(allowed.ncols()):
                v = {m: allowed[i, j] for i, m in enumerate(marks) if allowed[i, j]}
                row = {
                    i: value
                    for i, (f, m) in enumerate(monomials)
                    if (value := prod(u.get(t, 0) for t in f) * v.get(m, 0))
                }
                if not append(row):
                    continue
                direct_curve(
                    frames.colour,
                    q,
                    linear_combination(tangent, u),
                    linear_combination(two, v),
                )
                recorded.append(
                    dict(
                        tangent=sorted(u.items()),
                        mark=sorted((m, str(vv)) for m, vv in v.items()),
                        evaluation=sorted((i, str(x)) for i, x in row.items()),
                    )
                )
        if len(independent) == threshold:
            break
    # Verify containment of the incoming image in the kernel by native rational matrices.
    ev = fmpq_mat(len(recorded), n)
    for i, r in enumerate(recorded):
        for j, x in r["evaluation"]:
            ev[i, j] = fmpq(x)
    assert not any((ev * incoming).entries())
    stacked = array(
        n,
        [[ev[i, j] for j in range(n)] for i in range(ev.nrows())]
        + [[outgoing[i, j] for j in range(n)] for i in range(outgoing.nrows())],
    )
    assert stacked.rank() == len(independent) <= threshold
    return dict(
        root_colours=q,
        short_axis=short,
        order=order,
        initial_dimension=n,
        incoming_rank=incoming.rank(),
        outgoing_rank=outgoing.rank(),
        quadratic_cohomology=base["surviving_dimension"],
        closed_arcs_tried=tried,
        independent_evaluations=recorded,
        symbols=monomials,
        global_source_symbol_upper_bound=threshold - len(independent),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--quartet", nargs=4, type=int, default=(8, 19, 27, 47))
    p.add_argument("--short-axis", type=int, choices=range(4), default=1)
    p.add_argument("--order", type=int, default=6)
    p.add_argument("--samples", type=int, default=4)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if not __debug__ or a.samples < 1 or not 0 <= a.order <= 8:
        p.error("enable assertions, positive samples, and an order through eight")
    result = calculate(tuple(a.quartet), a.short_axis, a.order, a.samples)
    a.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        result["quadratic_cohomology"],
        len(result["independent_evaluations"]),
        result["global_source_symbol_upper_bound"],
    )
