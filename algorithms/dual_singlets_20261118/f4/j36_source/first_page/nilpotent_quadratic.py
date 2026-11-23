"""Test the quadratic commuting equations on a weight of the transverse F4 complex."""

import argparse
from collections import defaultdict
from functools import cache
from itertools import combinations_with_replacement
import json
from pathlib import Path

from flint import fmpq, fmpq_mat

from nilpotent_tangent import Colour, add, encode, simple_coordinates, weight_blocks


def pivots(matrix):
    reduced, rank = matrix.rref()
    return reduced, [
        next(j for j in range(reduced.ncols()) if reduced[i, j]) for i in range(rank)
    ]


def array(height, columns):
    result = fmpq_mat(height, len(columns))
    for j, col in enumerate(columns):
        for i, value in enumerate(col):
            if value:
                result[i, j] = value
    return result


def add_term(values, key, term):
    if term:
        values[key] = term + values.get(key, 0)
        if not values[key]:
            del values[key]


class Frames:
    def __init__(self, colour, roots, q):
        self.colour = colour
        self.forms = [[] for _ in range(5)]
        self.charts = {}
        self.splittings = {}
        for weight, raw, differentials in weight_blocks(colour, roots, q):
            ds = [fmpq_mat(m) for m in differentials]
            for degree in (1, 2, 3):
                basis = raw[degree]
                n = len(basis)
                if not n:
                    continue
                before = ds[degree - 1]
                _, incoming = pivots(before)
                red, columns = pivots(ds[degree])
                kernel = []
                for free in sorted(set(range(n)) - set(columns)):
                    v = [fmpq(0)] * n
                    v[free] = fmpq(1)
                    for i, j in enumerate(columns):
                        v[j] = -red[i, free]
                    kernel.append(v)
                images = [[before[i, j] for i in range(n)] for j in incoming]
                possible = images + kernel
                _, selected = pivots(array(n, possible))
                assert selected[: len(images)] == list(range(len(images))) and len(
                    selected
                ) == len(kernel)
                vectors = [possible[j] for j in selected]
                if degree == 2:
                    complements = [[int(k == j) for k in range(n)] for j in columns]
                    frame = array(n, vectors + complements)
                    assert frame.ncols() == n
                    self.splittings[weight] = (
                        basis,
                        frame.inv(),
                        raw[1],
                        incoming,
                        len(vectors) - len(images),
                        len(self.forms[2]),
                    )
                if not vectors:
                    continue
                full = array(n, vectors)
                _, rows = pivots(full.transpose())
                inverse = fmpq_mat([[v[i] for v in vectors] for i in rows]).inv()
                forms = self.forms[degree]
                labels = []
                for v in vectors[len(images) :]:
                    labels.append(len(forms))
                    forms.append((weight, {key: x for key, x in zip(basis, v) if x}))
                self.charts[degree, weight] = basis, rows, inverse, full, labels

    def bracket(self, i, degree, j):
        w, left = self.forms[1][i]
        v, right = self.forms[degree][j]
        tensor = {}
        for (axes, a), x in left.items():
            for (other, b), y in right.items():
                if axes[0] in other:
                    continue
                sign = (-1) ** sum(axes[0] > t for t in other)
                key = tuple(sorted((*axes, *other)))
                for k, z in self.colour.br[a][b].items():
                    add_term(tensor, (key, k), sign * x * y * z)
        if not tensor:
            return {}
        return self.project(degree + 1, add(w, v), tensor)

    def project(self, degree, weight, tensor):
        chart = self.charts[degree, weight]
        basis, rows, inverse, full, labels = chart
        target = array(len(basis), [[tensor.get(e, 0) for e in basis]])
        coordinates = inverse * array(len(rows), [[target[r, 0] for r in rows]])
        assert full * coordinates == target
        return {
            label: coordinates[i, 0]
            for i, label in enumerate(labels, full.ncols() - len(labels))
            if coordinates[i, 0]
        }


def symbol_monomials(frames, target, order):
    tangent = [encode(w) for w, _ in frames.forms[1]]
    marks = [encode(w) for w, _ in frames.forms[2]]
    max_digit = max(abs(x) for w, _ in frames.forms[1] + frames.forms[2] for x in w)
    assert (order + 1) * max_digit + max(map(abs, target)) < 1024
    half = order // 2
    left = defaultdict(list)
    for factors in combinations_with_replacement(range(len(tangent)), half):
        left[sum(tangent[i] for i in factors)].append(factors)
    wanted = -encode(target)
    result = []
    for right in combinations_with_replacement(range(len(tangent)), order - half):
        residue = wanted - sum(tangent[i] for i in right)
        for m, w in enumerate(marks):
            for factors in left.get(residue - w, ()):
                if not factors or not right or factors[-1] <= right[0]:
                    result.append((factors + right, m))
    return sorted(result)


def matrix_receipt(columns):
    keys = sorted({k for col in columns for k in col})
    matrix = array(len(keys), [[col.get(k, 0) for k in keys] for col in columns])
    _, picked = pivots(matrix)
    selected = array(len(keys), [[columns[j].get(k, 0) for k in keys] for j in picked])
    _, rows = pivots(selected.transpose())
    minor = fmpq_mat([[selected[i, j] for j in range(len(picked))] for i in rows])
    assert minor.det()
    return matrix, dict(
        rows=len(keys),
        columns=len(columns),
        rank=len(picked),
        selected_rows=[keys[i] for i in rows],
        selected_columns=picked,
        determinant=str(minor.det()),
    )


def calculate(q, short_axis, order, return_frames=False):
    c = Colour()
    roots = simple_coordinates(c)
    assert len(set(q)) == 4 and all(not c.br[i][j] for i in q for j in q)
    target = tuple(
        4 * s - 2 * r for s, r in zip(add(*(roots[i] for i in q)), roots[q[short_axis]])
    )
    frames = Frames(c, roots, q)
    columns = symbol_monomials(frames, target, order)
    coordinates = {key: i for i, key in enumerate(columns)}
    tangents, marks = len(frames.forms[1]), len(frames.forms[2])
    quadrics = [dict() for _ in range(marks)]
    reverse = {}
    for i, j in combinations_with_replacement(range(tangents), 2):
        values = frames.bracket(i, 1, j)
        reverse[i, j] = values
        for m, value in values.items():
            quadrics[m][i, j] = value / (1 + int(i == j))

    @cache
    def third(i, m):
        return frames.bracket(i, 2, m)

    outgoing = []
    pairs, triples = set(), set()
    for factors, m in columns:
        outgoing.append(
            {tuple(sorted((*factors, *edge))): v for edge, v in quadrics[m].items()}
        )
        for pos, i in enumerate(factors):
            spectators = factors[:pos] + factors[pos + 1 :]
            for z in third(i, m):
                triples.add((spectators, z))
            for other, j in enumerate(spectators):
                rest = spectators[:other] + spectators[other + 1 :]
                for n in reverse[min(i, j), max(i, j)]:
                    if n != m:
                        pairs.add((rest, min(m, n), max(m, n)))
    boundaries = []
    for factors, m, n in sorted(pairs):
        out = {}
        for active, mark, sign in ((m, n, 1), (n, m, -1)):
            for (i, j), v in quadrics[active].items():
                add_term(
                    out, coordinates[tuple(sorted((*factors, i, j))), mark], sign * v
                )
        boundaries.append(out)
    for factors, z in sorted(triples):
        out = {}
        # Project every compatible term in the three-form bracket.
        wz = frames.forms[3][z][0]
        for i, (wi, _) in enumerate(frames.forms[1]):
            for m, (wm, _) in enumerate(frames.forms[2]):
                if add(wi, wm) == wz:
                    value = third(i, m).get(z)
                    if value:
                        out[coordinates[tuple(sorted((*factors, i))), m]] = value
        boundaries.append(out)
    outgoing_matrix, out = matrix_receipt(outgoing)
    _, inc = matrix_receipt(boundaries)
    incoming_matrix = array(
        len(columns), [[b.get(i, 0) for i in range(len(columns))] for b in boundaries]
    )
    assert not any((outgoing_matrix * incoming_matrix).entries())
    result = dict(
        root_colours=q,
        short_axis=short_axis,
        order=order,
        cohomology_dimensions=[len(s) for s in frames.forms[1:4]],
        basis="rational RREF in native lexicographic weight, Taylor axes and colour order",
        source_symbols=columns,
        quadratic_outgoing=out,
        quadratic_incoming=inc,
        two_mark_parents=len(pairs),
        three_jet_parents=len(triples),
        surviving_dimension=len(columns) - out["rank"] - inc["rank"],
    )
    if return_frames:
        return result, frames, columns, incoming_matrix, outgoing_matrix
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quartet", nargs=4, type=int, default=(8, 19, 27, 47))
    parser.add_argument("--short-axis", type=int, choices=range(4), default=1)
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--out", type=Path, required=True)
    a = parser.parse_args()
    if not __debug__ or not 0 <= a.order <= 8:
        parser.error(
            "enable assertions and use a transverse degree from zero through eight"
        )
    data = calculate(tuple(a.quartet), a.short_axis, a.order)
    a.out.write_text(json.dumps(data, indent=2) + "\n")
    print(
        len(data["source_symbols"]),
        data["quadratic_incoming"]["rank"],
        data["quadratic_outgoing"]["rank"],
        data["surviving_dimension"],
    )
