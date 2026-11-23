"""Find exact polynomial commuting families and transported marks on the nilpotent F4 slices."""

import argparse
from collections import defaultdict
from itertools import combinations
import json
from math import prod
from pathlib import Path
from random import Random

from flint import fmpq, fmpq_mat
from nilpotent_arcs import linear_combination
from nilpotent_cubic import kernel, raw_bracket
from nilpotent_quadratic import add_term, array, calculate as quadratic, pivots
from nilpotent_tangent import simple_coordinates, weight_blocks


class Differential:
    def __init__(self, colour, roots, q):
        self.section = {}
        for _, spaces, maps in weight_blocks(colour, roots, q):
            for degree in (2, 3):
                d = fmpq_mat(maps[degree - 1])
                _, cols = pivots(d)
                basis = spaces[degree]
                image = array(
                    len(basis), [[d[i, j] for i in range(d.nrows())] for j in cols]
                )
                _, rows = pivots(image.transpose())
                inv = fmpq_mat(
                    [[image[i, j] for j in range(len(cols))] for i in rows]
                ).inv()
                for i, key in enumerate(basis):
                    h = {}
                    residual = {key: fmpq(1)}
                    if i in rows:
                        slot = rows.index(i)
                        for j, col in enumerate(cols):
                            x = inv[j, slot]
                            if x:
                                h[spaces[degree - 1][col]] = x
                                for k, target in enumerate(basis):
                                    add_term(residual, target, -d[k, col] * x)
                    self.section[key] = h, residual

    def reduce(self, form):
        preimage = {}
        residual = {}
        for key, value in form.items():
            h, r = self.section[key]
            for k, x in h.items():
                add_term(preimage, k, -value * x)
            for k, x in r.items():
                add_term(residual, k, value * x)
        return preimage, residual


def residual_at(c, curve, mark, t, q):
    # Test the original equations, including components outside every chosen cohomology frame.
    x = {((i,), k): 1 for i, k in enumerate(q)}
    y = {}
    for n, v in enumerate(curve[1:], 1):
        for k, a in v.items():
            add_term(x, k, t**n * a)
    for n, v in enumerate(mark):
        for k, a in v.items():
            add_term(y, k, t**n * a)
    return raw_bracket(c, x, x), raw_bracket(c, x, y)


def curve_for(c, differential, first, limit):
    # curve[n] is the coefficient of t^n beyond the fixed commuting origin.
    curve = [{}, first]
    for n in range(2, 2 * limit + 1):
        residual = {}
        for i in range(max(1, n - len(curve) + 1), min(n // 2, len(curve) - 1) + 1):
            j = n - i
            if j >= len(curve):
                continue
            for key, v in raw_bracket(c, curve[i], curve[j]).items():
                add_term(residual, key, v / (1 + int(i == j)))
        if n <= limit:
            correction, rem = differential.reduce(residual)
            if rem:
                return None
            curve.append(correction)
        elif residual:
            return None
    while not curve[-1]:
        curve.pop()
    return curve


def mark_columns(c, differential, curve, marks, limit, forms):
    columns = []
    extensions = []
    for m in marks:
        series = [forms[m][1]]
        col = {}
        for n in range(1, len(curve) + limit):
            residual = {}
            for i in range(1, min(n, len(curve) - 1) + 1):
                if n - i >= len(series):
                    continue
                for key, v in raw_bracket(c, curve[i], series[n - i]).items():
                    add_term(residual, key, v)
            if n <= limit:
                correction, residual = differential.reduce(residual)
                series.append(correction)
            for key, v in residual.items():
                col[n, key] = v
        columns.append(col)
        extensions.append(series)
    keys = sorted(set().union(*(c.keys() for c in columns)))
    system = array(len(keys), [[col.get(k, 0) for k in keys] for col in columns])
    return kernel(system), extensions


def calculate(q, short, order, previous, max_degree=3, samples=3):
    base, f, mons, inc, out = quadratic(q, short, order, True)
    assert (
        tuple(previous["root_colours"]) == q
        and previous["short_axis"] == short
        and previous["order"] == order
    )
    n = len(mons)
    threshold = n - inc.rank()
    # Extend the rows already certified on straight families.
    rows = []

    def append(row):
        row = {k: fmpq(x) for k, x in row.items() if x}
        for p, v in rows:
            a = row.get(p, 0)
            if a:
                for k, b in v.items():
                    add_term(row, k, -a * b)
        if not row:
            return False
        i = min(row)
        a = row[i]
        rows.append((i, {k: b / a for k, b in row.items()}))
        return True

    for i in range(out.nrows()):
        append({j: out[i, j] for j in range(n) if out[i, j]})
    for e in previous["independent_evaluations"]:
        append(dict(e["evaluation"]))
    assert n - len(rows) - inc.rank() == previous["global_source_symbol_upper_bound"]
    before = len(rows)
    result = dict(
        root_colours=q,
        short_axis=short,
        order=order,
        prior_bound=threshold - len(rows),
        max_polynomial_degree=max_degree,
        samples=samples,
        additional_evaluations=[],
    )
    if len(rows) == threshold:
        result["global_source_symbol_upper_bound"] = 0
        return result
    differential = Differential(f.colour, simple_coordinates(f.colour), q)
    by_mark = defaultdict(set)
    for factors, m in mons:
        by_mark[m].add(tuple(sorted(set(factors))))
    minimal = set().union(*by_mark.values())
    extra = {
        tuple(sorted(set(a) | set(b)))
        for v in by_mark.values()
        for a, b in combinations(v, 2)
    }
    supports = sorted(minimal, key=lambda s: (len(s), s)) + sorted(
        extra - minimal, key=lambda s: (len(s), s)
    )
    rng = Random(144363)
    tested = set()
    for axes in supports:
        if len(rows) == threshold:
            break
        for trial in range(samples):
            u = {i: (1 if trial == 0 else rng.choice((-2, -1, 1, 2))) for i in axes}
            key = tuple(sorted(u.items()))
            if key in tested:
                continue
            tested.add(key)
            curve = curve_for(
                f.colour, differential, linear_combination(f.forms[1], u), max_degree
            )
            if curve is None:
                continue
            # Components with zero leading evaluation may be required by the Bianchi equations.
            allowed, extensions = mark_columns(
                f.colour,
                differential,
                curve,
                range(len(f.forms[2])),
                max_degree,
                f.forms[2],
            )
            for j in range(allowed.ncols()):
                v = {m: allowed[m, j] for m in range(allowed.nrows()) if allowed[m, j]}
                row = {
                    k: b
                    for k, (factors, m) in enumerate(mons)
                    if (b := prod(u.get(s, 0) for s in factors) * v.get(m, 0))
                }
                if not append(row):
                    continue
                ys = [{} for _ in range(max_degree + 1)]
                for m, b in v.items():
                    for s, coef in enumerate(extensions[m]):
                        for k, x in coef.items():
                            add_term(ys[s], k, b * x)
                # The commuting polynomial has degree at most 2d; the marked equation also does.
                # Evaluating at 2d+1 different integers proves both identities exactly.
                for t in range(2 * max_degree + 1):
                    assert residual_at(f.colour, curve, ys, t, q) == ({}, {})
                result["additional_evaluations"].append(
                    dict(
                        tangent=sorted(u.items()),
                        mark=sorted((m, str(a)) for m, a in v.items()),
                        scalar_degree=len(curve) - 1,
                        mark_degree=max(i for i, coef in enumerate(ys) if coef),
                        evaluation=sorted((k, str(a)) for k, a in row.items()),
                    )
                )
            if len(rows) == threshold:
                break
    assert len(result["additional_evaluations"]) == len(rows) - before
    for rec in result["additional_evaluations"]:
        row = dict(rec["evaluation"])
        assert not any(
            (fmpq_mat([[fmpq(row.get(k, 0)) for k in range(n)]]) * inc).entries()
        )
    result["global_source_symbol_upper_bound"] = threshold - len(rows)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("previous", type=Path)
    parser.add_argument("--degree", type=int, default=3)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--out", required=True, type=Path)
    a = parser.parse_args()
    if not __debug__ or min(a.degree, a.samples) < 1:
        parser.error("enable assertions and positive degree and samples")
    prev = json.loads(a.previous.read_bytes())
    res = calculate(
        tuple(prev["root_colours"]),
        prev["short_axis"],
        prev["order"],
        prev,
        a.degree,
        a.samples,
    )
    a.out.write_text(json.dumps(res, indent=2) + "\n")
    print(res["prior_bound"], res["global_source_symbol_upper_bound"])
