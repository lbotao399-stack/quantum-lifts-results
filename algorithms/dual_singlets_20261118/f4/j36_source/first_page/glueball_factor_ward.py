"""Solve the shortened F4 terminal Ward equations at its highest Taylor weight."""

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, product
import json
from pathlib import Path

from flint import fmpq_mat, nmod_mat

from colour import Colour, Q, add, compact, items, prepend, wedge
from dotted_filtration import pure_current
from intermediate_rank import tensor_multiplicity
from source_symmetry import (
    centralizer,
    character,
    multiplicity,
    partitions,
    singlet_trace,
)


def negative_sum(*weights):
    return tuple(-sum(terms) for terms in zip(*weights))


def product_types():
    full = (4, 4, 3, 3)
    types = []
    for partition in partitions(12):
        mu = partition[::-1]
        mult = sum(
            Fraction(
                character(mu, a)
                * character((1, 1), b)
                * character(full, tuple(sorted(a + b))),
                centralizer(a) * centralizer(b),
            )
            for a in partitions(12)
            for b in partitions(2)
        )
        if not mult:
            continue
        assert mult == 1
        # Only a full third row of three boxes permits the first Berezinian shift.
        after = None
        if len(mu) >= 3 and mu[2] == 3:
            after = (mu[0] + 1, mu[1] + 1, *mu[3:])
        types.append(
            dict(
                shortened_type=mu,
                first_target=after,
                source_cochains=multiplicity(mu, 11, singlet_trace),
                terminal_target_cochains=multiplicity(after, 11, singlet_trace)
                if after
                else 0,
            )
        )
    assert {r["shortened_type"] for r in types} == {
        (4, 4, 2, 2),
        (4, 3, 3, 2),
        (3, 3, 3, 3),
    }
    return types


class Anchors:
    def __init__(self, colour, highest):
        self.c = colour
        self.h = highest
        self.low = colour.low[4, 0]
        self.hweight = colour.weight(highest)
        self.lweight = colour.weight(self.low)
        target_weight = negative_sum(self.hweight, self.lweight)
        self.scalar_pairs = [
            wedge(pair)[0]
            for pair in combinations(range(52), 2)
            if colour.weight(wedge(pair)[0]) == target_weight
        ]
        self.index = {pair: i for i, pair in enumerate(self.scalar_pairs)}
        self.pairs = defaultdict(list)
        for i, j in product(range(52), repeat=2):
            self.pairs[
                tuple(x + y for x, y in zip(colour.weights[i], colour.weights[j]))
            ].append((i, j))

    def at_low(self, t):
        return {self.index[t]: Q(1)}

    @cache
    def evaluate(self, v, t):
        c = self.c
        w = c.weight(v)
        if c.weight(t) != negative_sum(w, self.hweight):
            return {}
        if w == self.lweight:
            # Uniqueness of this extremal wedge in the full exterior power is checked below.
            if v != self.low:
                raise AssertionError("a second wedge has the extremal weight")
            return self.at_low(t)
        shifted = tuple(x - y for x, y in zip(w, c.rho2))
        den = 18 * self.low.bit_count() - c.scalar_product(w, shifted)
        if den <= 0:
            return {}
        row = {}
        for e, f, g in c.roots:
            for vv, x in c.act(f, v).items():
                for tt, y in c.act(e, t).items():
                    for j, z in self.evaluate(vv, tt).items():
                        add(row, j, g * x * y * z)
        return {j: -2 * x / den for j, x in row.items() if x}

    def mixed(self, three, scalar, marked):
        row = {}
        for pair, coefficient in self.c.delta[marked].items():
            x, y = items(pair)
            for i, j, s in ((x, y, 1), (y, x, -1)):
                v, sv = prepend(i, three)
                t, st = prepend(j, 1 << scalar)
                if sv and st:
                    for k, z in self.evaluate(v, t).items():
                        add(row, k, sv * st * s * coefficient * z)
        return compact(row)

    def constraints(self, depth_limit):
        c = self.c
        # The scalar-scalar bracket and the stabilizer of the two extremal wedges.
        for i in range(52):
            row = {self.index[t]: x for t, x in c.delta[i].items() if t in self.index}
            if row:
                yield ("scalar", i), row
            if c.act(i, self.h) or c.act(i, self.low):
                continue
            for a, b in self.pairs[
                negative_sum(self.hweight, self.lweight, c.weights[i])
            ]:
                if a >= b:
                    continue
                row = {
                    self.index[t]: x
                    for t, x in c.act(i, wedge((a, b))[0]).items()
                    if t in self.index
                }
                if row:
                    yield ("stabilizer", i, a, b), row
        active = {self.low ^ (1 << i) for i in items(self.low)}
        seen = set()
        for _ in range(depth_limit + 1):
            new = set()
            for three in sorted(active - seen):
                for marked, scalar in self.pairs[
                    negative_sum(self.hweight, c.weight(three))
                ]:
                    row = self.mixed(three, scalar, marked)
                    if row:
                        yield ("mixed", items(three), scalar, marked), row
                for e in c.data["raising"]:
                    new.update(c.act(e, three))
            seen.update(active)
            active = new - seen
        active = {self.low}
        seen = set()
        for _ in range(depth_limit + 1):
            new = set()
            for v in sorted(active - seen):
                weight = negative_sum(self.hweight, c.weight(v))
                for r, w in enumerate(c.weights):
                    if w != weight:
                        continue
                    row = {}
                    for t, x in c.delta[r].items():
                        for k, z in self.evaluate(v, t).items():
                            add(row, k, x * z)
                    if compact(row):
                        yield ("scalar_at", items(v), r), compact(row)
                for e in c.data["raising"]:
                    raised = c.act(e, v)
                    new.update(raised)
                    total = negative_sum(self.hweight, c.weight(v), c.weights[e])
                    for a, b in self.pairs[total]:
                        if a >= b:
                            continue
                        t = wedge((a, b))[0]
                        row = {}
                        for vv, x in raised.items():
                            for k, z in self.evaluate(vv, t).items():
                                add(row, k, x * z)
                        for tt, x in c.act(e, t).items():
                            for k, z in self.evaluate(v, tt).items():
                                add(row, k, x * z)
                        if compact(row):
                            yield ("raising", items(v), e, a, b), compact(row)
            seen.update(active)
            active = new - seen

    def solve(self, depth_limit):
        return solve_constraints(self.constraints(depth_limit), len(self.index))


def solve_constraints(equations, width):
    pivots = {}
    tested = 0
    witnesses = []
    matrix = []
    for probe, row in equations:
        tested += 1
        reduced = row.copy()
        for pivot, old in sorted(pivots.items()):
            coeff = reduced.get(pivot, 0)
            if coeff:
                for k, v in old.items():
                    add(reduced, k, -coeff * v)
        reduced = compact(reduced)
        if reduced:
            p = min(reduced)
            pivots[p] = {k: Q(v) / reduced[p] for k, v in reduced.items()}
            witnesses.append(probe)
            matrix.append([row.get(i, 0) for i in range(width)])
            if len(pivots) == width:
                break
    columns = sorted(pivots)
    det = fmpq_mat([[row[j] for j in columns] for row in matrix]).det()
    assert det
    modular = {}
    for prime in (1009, 1013):
        residual = [
            [
                int(Q(v).numerator) * pow(int(Q(v).denominator), -1, prime) % prime
                for v in row
            ]
            for row in matrix
        ]
        minor = nmod_mat([[row[j] for j in columns] for row in residual], prime)
        modular[str(prime)] = int(minor.det())
        assert (
            int(det.numerator) - modular[str(prime)] * int(det.denominator)
        ) % prime == 0
    return dict(
        unknowns=width,
        rank=len(pivots),
        nonempty_equations_visited=tested,
        probes=witnesses,
        pivot_columns=columns,
        rational_rows=[[str(x) for x in row] for row in matrix],
        determinant=str(det),
        minor_residues=modular,
    )


class BothLegs:
    def __init__(self, c, solvers):
        self.c = c
        self.solvers = {s.h: s for s in solvers}
        self.offset = {
            s.h: sum(len(v.index) for v in solvers[:i]) for i, s in enumerate(solvers)
        }
        self.high_weights = {c.weight(u): u for u in self.solvers}
        self.low = c.low[4, 0]
        self.width = sum(len(s.index) for s in solvers)

    @cache
    def evaluate(self, u, v, t):
        c = self.c
        w = c.weight(u)
        if c.weight(t) != negative_sum(w, c.weight(v)):
            return {}
        if w in self.high_weights:
            if u != self.high_weights[w]:
                raise AssertionError("ambiguous extremal wedge")
            offset = self.offset[u]
            return {j + offset: x for j, x in self.solvers[u].evaluate(v, t).items()}
        den = 18 * u.bit_count() - c.scalar_product(w, tuple(x + y for x, y in zip(w, c.rho2)))
        if den <= 0:
            return {}
        row = {}
        for e, f, g in c.roots:
            for uu, x in c.act(e, u).items():
                for vv, y in c.act(f, v).items():
                    for k, z in self.evaluate(uu, vv, t).items():
                        add(row, k, g * x * y * z)
                for tt, y in c.act(f, t).items():
                    for k, z in self.evaluate(uu, v, tt).items():
                        add(row, k, g * x * y * z)
        return {k: -2 * x / den for k, x in row.items() if x}

    def constraints(self, depth, lower):
        c = self.c
        for solver, block in zip(self.solvers.values(), lower):
            for probe, row in zip(block["probes"], block["rational_rows"]):
                yield (
                    ("lower", items(solver.h), probe),
                    {
                        i + self.offset[solver.h]: Q(x)
                        for i, x in enumerate(row)
                        if Q(x)
                    },
                )
        pairs = next(iter(self.solvers.values())).pairs
        active = {h ^ (1 << i) for h in self.solvers for i in items(h)}
        seen = set()
        for _ in range(depth + 1):
            new = set()
            for four in sorted(active - seen):
                for marked, scalar in pairs[
                    negative_sum(c.weight(four), c.weight(self.low))
                ]:
                    row = {}
                    for pair, coefficient in c.delta[marked].items():
                        x, y = items(pair)
                        for i, j, s in ((x, y, 1), (y, x, -1)):
                            u, su = prepend(i, four)
                            t, st = prepend(j, 1 << scalar)
                            if su and st:
                                for k, z in self.evaluate(u, self.low, t).items():
                                    add(row, k, su * st * s * coefficient * z)
                    if compact(row):
                        yield (
                            ("positive_mixed", items(four), scalar, marked),
                            compact(row),
                        )
                for f in c.data["lowering"]:
                    new.update(c.act(f, four))
            seen.update(active)
            active = new - seen
        for h in self.solvers:
            for hi in items(h):
                four = h ^ (1 << hi)
                for li in items(self.low):
                    three = self.low ^ (1 << li)
                    for r, w in enumerate(c.weights):
                        for a, b in pairs[
                            negative_sum(c.weight(four), c.weight(three), w)
                        ]:
                            if a >= b:
                                continue
                            t = wedge((a, b))[0]
                            row = {}
                            for pair, coefficient in c.delta[r].items():
                                x, y = items(pair)
                                for i, j, phase in ((x, y, 1), (y, x, -1)):
                                    u, su = prepend(i, four)
                                    v, sv = prepend(j, three)
                                    if su and sv:
                                        for k, z in self.evaluate(u, v, t).items():
                                            add(
                                                row,
                                                k,
                                                su * sv * phase * coefficient * z,
                                            )
                            if compact(row):
                                yield (
                                    ("ghost_pair", items(four), items(three), r, a, b),
                                    compact(row),
                                )


def calculate(depth):
    c = Colour()
    lowest = c.low[4, 0]
    target = c.weight(lowest)
    # Count every wedge of weight target by two ordered halves.
    pairs = defaultdict(list)
    for a, b in combinations(range(52), 2):
        weight = tuple(x + y for x, y in zip(c.weights[a], c.weights[b]))
        pairs[weight].append((a, b))
    at_low = [
        (a, b, d, e)
        for w, s in list(pairs.items())
        for a, b in s
        for d, e in pairs.get(tuple(x - y for x, y in zip(target, w)), ())
        if b < d
    ]
    assert at_low == [items(lowest)]
    for h in c.high[5]:
        w = c.weight(h)
        at_high = []
        for d, e, f in combinations(range(52), 3):
            triple_weight = tuple(
                x + y + z for x, y, z in zip(c.weights[d], c.weights[e], c.weights[f])
            )
            needed = tuple(x - y for x, y in zip(w, triple_weight))
            at_high.extend((a, b, d, e, f) for a, b in pairs.get(needed, ()) if b < d)
        assert at_high == [items(h)]
    (lows,) = pure_current(4, 0)
    (scalar,) = pure_current(2, 0)
    blocks = []
    solvers = []
    for labels, high in zip(sorted(pure_current(5, 0), reverse=True), c.high[5]):
        solver = Anchors(c, high)
        solvers.append(solver)
        record = solver.solve(depth)
        record.update(
            dynkin=labels,
            highest=items(high),
            lowest=items(lowest),
            scalar_pair_order=list(map(items, solver.scalar_pairs)),
            harmonic_tensor_multiplicity=tensor_multiplicity(labels, lows, scalar),
        )
        blocks.append(record)
    joint = BothLegs(c, solvers)
    equations = solve_constraints(joint.constraints(depth, blocks), joint.width)
    # Separate ordinary invariant-tensor conditions from the actual mixed brackets.
    conditions = []
    brackets = []
    for probe, row in zip(equations["probes"], equations["rational_rows"]):
        constraints = (
            conditions
            if probe[0] == "lower" and probe[2][0] in ("scalar", "stabilizer")
            else brackets
        )
        constraints.append(row)
    fixed = fmpq_mat(conditions)
    rr, rank = fixed.rref()
    pivots = [next(j for j in range(rr.ncols()) if rr[i, j]) for i in range(rank)]
    free = [j for j in range(rr.ncols()) if j not in pivots]
    change = fmpq_mat(rr.ncols(), len(free))
    for k, j in enumerate(free):
        change[j, k] = 1
        for i, p in enumerate(pivots):
            change[p, k] = -rr[i, j]
    assert fixed * change == fmpq_mat(fixed.nrows(), change.ncols())
    reduced = fmpq_mat(brackets) * change
    assert change.ncols() == sum(
        block["harmonic_tensor_multiplicity"] for block in blocks
    )
    assert reduced.nrows() == reduced.ncols() and reduced.det()
    control = dict(
        harmonic_dimension=change.ncols(),
        bracket_rows=reduced.tolist(),
        determinant=str(reduced.det()),
        basis=change.tolist(),
    )
    for key in ("bracket_rows", "basis"):
        control[key] = [list(map(str, row)) for row in control[key]]
    return dict(
        taylor_shape=[5, 4, 2],
        highest_physical_weight=[5, 4, 1, 1, 0],
        product_types=product_types(),
        negative_ghost_harmonic=lows,
        scalar_harmonic=scalar,
        blocks=blocks,
        joint=equations,
        harmonic_control=control,
        upper_bound=joint.width - equations["rank"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--out", type=Path, required=True)
    a = parser.parse_args()
    if not __debug__ or a.depth < 0:
        parser.error("enable assertions and use nonnegative depth")
    result = calculate(a.depth)
    a.out.write_text(json.dumps(result, indent=2) + "\n")
    print(
        [(b["unknowns"], b["rank"], b["determinant"]) for b in result["blocks"]],
        result["upper_bound"],
    )
