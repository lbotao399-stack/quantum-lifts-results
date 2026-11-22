"""Coinvariants of auxiliary sl(2) on the orthogonal determinant directions."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import sys
from flint import fmpq_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scalar_multipliers"))
from verify import path_types  # noqa: E402 — shared trace workspace above

WEIGHT = (0, 1, -1, 0)
RAISE = {0: (1, -2), 2: (0, 1)}
DUAL_RAISE = {0: (2, -1), 1: (0, 2)}
SHAPES = path_types(5, 5, 8, 4)


def canonical(raw, shape):
    atoms, plus, minus, extra = raw
    atoms = list(atoms)
    sign = 1
    for start, (closed, path) in enumerate(shape):
        if start and shape[start - 1] == (closed, path):
            continue
        stop = start + shape[start:].count((closed, path))
        seq = atoms[start:stop]
        if (path.count("m") + path.count("p")) % 2:
            if len(set(seq)) != len(seq):
                return 0, None
            sign *= (-1) ** sum(a > b for j, a in enumerate(seq) for b in seq[j + 1 :])
        atoms[start:stop] = sorted(seq)
    groups = []
    for seq in (plus, minus):
        if len(set(seq)) != len(seq):
            return 0, None
        sign *= (-1) ** sum(a > b for j, a in enumerate(seq) for b in seq[j + 1 :])
        groups.append(tuple(sorted(seq)))
    return sign, (tuple(atoms), *groups, tuple(sorted(extra)))


def directions(index, weight, shapes=SHAPES, extra_axes=2):
    shape, np, nm, _ = shapes[index]
    choices = []
    for closed, v in shape:
        axes = (3,) if closed or len(v) % 2 == 0 else range(3)
        choices.append(
            tuple(
                (b, ms) for b in axes for ms in product(range(3), repeat=v.count("s"))
            )
        )
    for atoms in product(*choices):
        if any(
            a > b
            for i, (a, b) in enumerate(zip(atoms, atoms[1:]))
            if shape[i] == shape[i + 1]
        ):
            continue
        q = sum(WEIGHT[b] - sum(WEIGHT[k] for k in m) for b, m in atoms)
        for plus, minus in product(
            combinations(range(3), np), combinations(range(3), nm)
        ):
            h = q + sum(WEIGHT[k] for k in plus + minus)
            for x in combinations_with_replacement(range(3), extra_axes):
                if h + sum(WEIGHT[k] for k in x) != weight:
                    continue
                raw = atoms, plus, minus, x
                sign, v = canonical(raw, shape)
                if sign:
                    yield v


def replace(seq, i, new):
    return seq[:i] + (new,) + seq[i + 1 :]


def raised(raw, shape):
    atoms, plus, minus, x = raw
    out = Counter()

    def add(v, c):
        sign, w = canonical(v, shape)
        if sign:
            out[w] += c * sign

    for i, (b, axes) in enumerate(atoms):
        if b in RAISE:
            new, n = RAISE[b]
            add((replace(atoms, i, (new, axes)), plus, minus, x), n)
        for j, a in enumerate(axes):
            if a in DUAL_RAISE:
                new, n = DUAL_RAISE[a]
                add((replace(atoms, i, (b, replace(axes, j, new))), plus, minus, x), n)
    for group in range(1, 4):
        for i, a in enumerate(raw[group]):
            if a in RAISE:
                new, n = RAISE[a]
                add(replace(raw, group, replace(raw[group], i, new)), n)
    return out


@lru_cache(None)
def quotient(i, shapes=SHAPES, extra_axes=2):
    shape = shapes[i][0]
    zero = list(directions(i, 0, shapes, extra_axes))
    lower = list(directions(i, -1, shapes, extra_axes))
    assert len(zero) == len(set(zero)) and len(lower) == len(set(lower))
    positions = {v: j for j, v in enumerate(zero)}
    m = fmpq_mat(len(lower), len(zero))
    for j, source in enumerate(lower):
        for v, n in raised(source, shape).items():
            m[j, positions[v]] = n
    red, rank = m.rref()
    assert rank == len(lower)
    pivots = {next(k for k in range(len(zero)) if red[j, k]) for j in range(rank)}
    return [j for j in range(len(zero)) if j not in pivots], zero
