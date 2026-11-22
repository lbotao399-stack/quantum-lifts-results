"""Rebuild the homogeneous annihilator of the first correction to the F4 sixth glueball."""
import argparse
from fractions import Fraction
from functools import cache
from itertools import combinations
import json
from random import Random

import numpy as np

from colour import items, prepend, wedge
from dotted_filtration import GRAM2, SIMPLE, casimir2, difference, dot, pure_current, root_coordinates
from f4_singlet_characters import reflect, roots
from source_symmetry import centralizer, character, multiplicity, partitions, singlet_trace
from ward import Descent


RHO = (1, 1, 1, 1)
ZERO = (0, 0, 0, 0)


@cache
def dominant(w):
    while min(w) < 0:
        i = next(i for i, v in enumerate(w) if v < 0)
        w = reflect(w, (i, SIMPLE[i]))
    return w


@cache
def orbit(w):
    entries = {w: 1}
    pending = [w]
    while pending:
        current = pending.pop()
        for mirror in enumerate(SIMPLE):
            target = reflect(current, mirror)
            if target not in entries:
                entries[target] = -entries[current]
                pending.append(target)
            else:
                assert entries[target] == -entries[current]
    assert len(entries) == 1152
    return entries


@cache
def positive_roots():
    positive = [(r, root_coordinates(r)) for r in roots("F4") if min(root_coordinates(r)) >= 0]
    assert len(positive) == 24
    return positive


@cache
def weight_size(high, w):
    top = dominant(w)
    if w != top:
        return weight_size(high, top)
    room = root_coordinates(difference(high, w))
    if min(room) < 0:
        return 0
    if w == high:
        return 1
    numerator = 0
    for root, simple in positive_roots():
        limit = min(k // v for k, v in zip(room, simple) if v)
        dual = tuple(dot(col, root) for col in GRAM2)
        for j in range(1, limit + 1):
            higher = tuple(x + j * y for x, y in zip(w, root))
            numerator += 2 * dot(higher, dual) * weight_size(high, higher)
    denominator = casimir2(high) - casimir2(w)
    assert denominator > 0 and numerator % denominator == 0
    return numerator // denominator


def shifted(w):
    return tuple(i + 1 for i in w)


def tensor_multiplicity(u, v, k):
    return sum(sign * weight_size(k, difference(shifted(v), w))
               for w, sign in orbit(shifted(u)).items())


def exchange_trace(u, k):
    total = 0
    for w, sign in orbit(RHO).items():
        delta = difference(shifted(k), w)
        if all(x % 2 == 0 for x in delta):
            total += sign * weight_size(u, tuple(x // 2 for x in delta))
    return total


def characters():
    k, = pure_current(3, 0)
    two = sorted(pure_current(5, 0), reverse=True)
    assert len(two) == 2
    adjoint = (1, 0, 0, 0)
    assert tensor_multiplicity(adjoint, adjoint, adjoint) == 1
    assert exchange_trace(adjoint, adjoint) == -1
    blocks = []
    for i, u in enumerate(two):
        for v in two[i:]:
            n = tensor_multiplicity(u, v, k)
            swap = exchange_trace(u, k) if u == v else None
            size = (n + swap) // 2 if u == v else n
            blocks.append(dict(first=u, second=v, multiplicity=n, exchange=swap, dimension=size))
    assert [r["dimension"] for r in blocks] == [18, 13, 9]
    diagrams = []
    for part in partitions(13):
        shape = part[::-1]
        amount = sum(Fraction(character((5,5), a) * character((3,), b)
                              * character(shape, tuple(sorted(a+b))), centralizer(a) * centralizer(b))
                     for a in partitions(10) for b in partitions(3))
        if amount:
            assert amount == 1
            diagrams.append(dict(shape=shape, dimension=multiplicity(shape, 13, singlet_trace)))
    assert sum(r["dimension"] for r in diagrams) == 633
    return dict(scalar_harmonic=k, ghost_harmonics=two, blocks=blocks, terminal_diagrams=diagrams)


class Columns(Descent):
    """All forty native seed functionals on one raw tree boundary at a time."""

    def __init__(self, prime):
        super().__init__()
        self.prime = prime
        self.zero = np.zeros(len(self.seeds), dtype=np.int64)
        self.anchors = {}
        for column, (i, j, seed) in enumerate(self.seeds):
            for mask, value in seed.items():
                self.anchors.setdefault((i, j, mask), self.zero.copy())[column] = self.residue(value)

    def residue(self, x):
        return int(x.numerator) * pow(int(x.denominator), -1, self.prime) % self.prime

    @cache
    def anchor(self, i, j, mask):
        if (i, j) == (0, 1):
            other, sign = self.transformed(mask)
            phase = self.residue(sign * self.phase[0] / self.phase[1])
            return phase * self.anchor(1, 0, other) % self.prime
        out = self.zero.copy()
        for key, v in self.dual(mask).items():
            out += self.residue(v) * self.anchors.get((i, j, key), self.zero)
        return out % self.prime

    @cache
    def descend(self, i, j, v, t):
        c = self.c
        w = c.weight(v)
        if any(x+y+z for x,y,z in zip(w, c.weight(t), self.high_weight[5][i])):
            return self.zero
        if w in self.low_weight[5]:
            if v == c.low[5,j]:
                return c.low_sign[5,j] * self.anchor(i, j, t) % self.prime
            return self.zero
        factor = self.denominator(5, w, False)
        if not factor:
            return self.zero
        out = self.zero.copy()
        for e, f, g in c.roots:
            for vv, x in c.act(f, v).items():
                for tt, y in c.act(e, t).items():
                    out += self.residue(g*x*y) * self.descend(i, j, vv, tt)
        return (out % self.prime) * self.residue(factor) % self.prime

    def row(self, probe):
        i, r, a, b = probe
        base = self.c.low[4, 0]
        active = wedge((a, b))[0]
        out = self.zero.copy()
        for pair, z in self.c.delta[r].items():
            x, y = items(pair)
            for first, second, factor in ((x, y, z), (y, x, -z)):
                v, sv = prepend(first, base)
                t, st = prepend(second, active)
                if sv and st:
                    coeff = self.residue(factor * sv * st)
                    out += coeff * (self.descend(i, 0, v, t) + self.descend(i, 1, v, t))
        return out % self.prime

    def probes(self):
        c = self.c
        pairs = {}
        for a, b in combinations(range(52), 2):
            pairs.setdefault(c.weight(wedge((a, b))[0]), []).append((a, b))
        rows = []
        for i, h in enumerate(c.high[5]):
            for r, w in enumerate(c.weights):
                need = tuple(-x-y-z for x,y,z in zip(c.weight(h), c.weight(c.low[4,0]), w))
                rows.extend((i, r, a, b) for a, b in pairs.get(need, ()))
        assert len(rows) == 763
        Random(4433633).shuffle(rows)
        return rows


def rank_check(prime, probes=None):
    columns = Columns(prime)
    basis = {}
    selected = []
    matrix = []
    for probe in probes if probes is not None else columns.probes():
        row = columns.row(probe)
        reduced = row.copy()
        for pivot in sorted(basis):
            reduced = (reduced - reduced[pivot] * basis[pivot]) % prime
        indices = np.flatnonzero(reduced)
        if len(indices):
            pivot = int(indices[0])
            basis[pivot] = reduced * pow(int(reduced[pivot]), -1, prime) % prime
            matrix.append(row.tolist())
            selected.append(probe)
            if len(basis) == 40:
                break
    from sympy import Matrix
    assert len(selected) == 40
    det = int(Matrix(matrix).det(method="domain-ge")) % prime
    assert det
    columns.anchor.cache_clear()
    columns.descend.cache_clear()
    return dict(prime=prime, rank=len(selected), determinant=det), selected


def calculate():
    dims = characters()
    first, probes = rank_check(1000003)
    second, check = rank_check(1000033, probes)
    assert check == probes
    return dict(characters=dims, probes=probes, minors=[first, second],
                intermediate_classical_dimension=0, projected_intermediate_classical_dimension=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--characters", action="store_true")
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    print(json.dumps(characters() if args.characters else calculate(), indent=2))
