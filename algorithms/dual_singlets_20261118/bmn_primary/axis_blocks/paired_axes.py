"""Bound a symplectic Taylor cohomology before computing the physical differential."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, product
import json
from math import factorial

from transverse import (PRIMARY, character, colour_data, low_character, partitions,
                        schur_table, signed_orbit, weight_tables)
from b2_schur_complex import centralizer, multiplicity, partner
from bmn_primary_index import G2_WEYL
from symplectic_axis_colours import singlet_bound


def axis_chars():
    records = json.loads((PRIMARY / "symplectic_655/symplectic_axis_bounds.json").read_bytes())["blocks"]
    weights = defaultdict(dict)
    for energy, (x, y), length, _, _, _, dimension in records:
        if dimension:
            for w in signed_orbit(x, y):
                weights[energy, length][w] = dimension
    return records, dict(weights)


AXIS_RECORDS, AXIS = axis_chars()

@cache
def axis(group):
    if group == "USp(4)":
        return AXIS
    result = {(0, 0): {(0, 0): 1}}
    for energy in range(1, 8):
        for (x, y, length), n in weight_tables((energy,), group)[1].items():
            result.setdefault((energy, length), {}).update({v: n for v in signed_orbit(x, y, group)})
    return result


@cache
def exterior(group, power):
    result = Counter()
    weights, _ = colour_data(group)
    for selection in combinations(weights, power):
        result[sum(x for x, _ in selection), sum(y for _, y in selection)] += 1
    return result


@cache
def factor_character(group, rows, factor):
    name, orbit, degree, length = factor
    if name == "low":
        return low_character(orbit, degree, length, group, rows)
    weights = axis(group)[degree, length] if name == "axis" else exterior(group, degree)
    return {(orbit * x, orbit * y): count for (x, y), count in weights.items()}


@cache
def tensor(group, rows, factors):
    if not factors:
        return {(0, 0): 1}
    result = Counter()
    for (a, b), n in tensor(group, rows, factors[:-1]).items():
        for (x, y), m in factor_character(group, rows, factors[-1]).items():
            result[a + x, b + y] += n * m
    return {w: n for w, n in result.items() if n}


def trivial_multiplicity(weights, group):
    if group == "G2":
        return sum(s * weights.get((1 - x, 1 - y), 0) for (x, y), s in G2_WEYL.items())
    return sum(sign * i * j * weights.get((2 - i * a, 1 - j * b), 0)
               for (a, b), sign in (((2, 1), 1), ((1, 2), -1))
               for i, j in product((-1, 1), repeat=2))

@cache
def singlets(group, rows, factors):
    return trivial_multiplicity(tensor(group, rows, factors), group)


def permutation(cycles):
    perm, masks = [], []
    for length in cycles:
        start = len(perm)
        perm.extend(range(start + 1, start + length))
        perm.append(start)
        masks.append(((1 << length) - 1) << start)
    return tuple(perm), masks


def submasks(mask):
    current = mask
    while current:
        yield current
        current = (current - 1) & mask
    yield 0


@cache
def mixed_patterns(a, b, mask, length, swap, perm):
    """Orbit-stable occupations of jets meeting at least two of the three blocks."""
    if length == 0:
        return {(): 1} if (a, b, mask) == (0, 0, 0) else {}
    if 2 * length > a + b + mask.bit_count() or length < 0:
        return {}
    jets = {(x, y, z) for x in range(a + 1) for y in range(b + 1) for z in submasks(mask)
            if (x > 0) + (y > 0) + (z > 0) >= 2}
    orbits = []
    while jets:
        jet = min(jets)
        dx = dy = low_mask = period = 0
        overlap = False
        while jet in jets:
            jets.remove(jet)
            x, y, z = jet
            dx += x
            dy += y
            overlap |= bool(z & low_mask)
            low_mask |= z
            period += 1
            moved = sum(1 << j for i, j in enumerate(perm) if z & (1 << i))
            jet = (y, x, moved) if swap else (x, y, moved)
        if not overlap and dx <= a and dy <= b and period <= length:
            orbits.append((dx, dy, low_mask, period))
    answer = Counter()

    def visit(start, left, right, remaining_mask, remaining, factors, sign):
        if remaining == 0:
            if (left, right, remaining_mask) == (0, 0, 0):
                answer[tuple(sorted(factors))] += sign
            return
        if left + right + remaining_mask.bit_count() < 2 * remaining:
            return
        for j in range(start, len(orbits)):
            dx, dy, z, period = orbits[j]
            if z & remaining_mask != z:
                continue
            limit = min(remaining // period, 1 if z else remaining,
                        left // dx if dx else remaining, right // dy if dy else remaining)
            for n in range(1, limit + 1):
                visit(j + 1, left - n * dx, right - n * dy, remaining_mask ^ z,
                      remaining - n * period, (*factors, ("mixed", period, n, 0)),
                      sign * (-1) ** ((period - 1) * n))

    visit(0, a, b, mask, length, (), 1)
    return {k: n for k, n in answer.items() if n}


@cache
def page_trace(high, cycles, length, swap, low_homology, group="USp(4)", rows=None):
    perm, cycle_masks = permutation(cycles)
    all_low = (1 << len(perm)) - 1
    pure = sorted((e, y) for e, y in axis(group) if e <= high)
    total = 0
    for (e, y), (f, t) in product(pure, repeat=2):
        if swap and (e, y) != (f, t):
            continue
        axes = []
        if swap:
            if e:
                axes.append(("axis", 2, e, y))
        else:
            axes = [("axis", 1, k, n) for k, n in ((e, y), (f, t)) if k]
        sign = (-1) ** y if swap else 1
        for selection in range(1 << len(cycles)):
            chosen = tuple(c for i, c in enumerate(cycles) if selection & (1 << i))
            mask = all_low ^ sum(m for i, m in enumerate(cycle_masks) if selection & (1 << i))
            for low_length in range(sum(chosen) + 1):
                degree_left = length - y - t - low_length
                if degree_left < 0 or 2 * degree_left > 2 * high - e - f + mask.bit_count():
                    continue
                if not low_character(chosen, low_length, low_homology, group, rows):
                    continue
                others = axes + ([("low", chosen, low_length, low_homology)] if chosen else [])
                for mixed, n in mixed_patterns(high - e, high - f, mask,
                                                degree_left, swap, perm).items():
                    total += sign * n * singlets(group, rows, tuple(sorted((*others, *mixed))))
    return total


def projection(high, tail, length, homology, group="USp(4)", rows=None):
    numerator = 0
    traces = []
    for cycles in partitions(sum(tail)):
        size = factorial(sum(tail)) // centralizer(cycles)
        specht = character(tail, cycles)
        values = [page_trace(high, cycles, length, swap, homology, group, rows) for swap in (False, True)]
        numerator += size * specht * (values[0] + (-1) ** high * values[1])
        traces.append([list(cycles), size, specht, *values])
    denominator = 2 * factorial(sum(tail))
    assert numerator >= 0 and numerator % denominator == 0
    return numerator // denominator, traces


def check_controls():
    original = singlet_bound(AXIS_RECORDS)
    for k in range(1, 5):
        symmetric = schur_table(k, True)[(k,)]
        for length in range(k + 1):
            actual = {(x, y): value for (x, y, n), value in symmetric.items() if n == length}
            expected = {w: n for w, n in AXIS.get((k, length), {}).items() if w[0] >= w[1] >= 0}
            assert actual == expected
    # Compare with a complete, independently projected four-slot calculation.
    degrees = []
    for length in range(5):
        p, _ = projection(1, (1, 1), length, True)
        # Alternation on both ordered pairs: s_(1,1) s_(1,1).
        components = ((2, 2), (2, 1, 1), (1, 1, 1, 1))
        cohom = sum(colour_singlets(schur_table(4, True)[s], length) for s in components)
        assert p >= cohom
        degrees.append([length, p, cohom])
        assert page_trace(1, (1, 1), length, False, False) == original((1, 1, 1, 1), length)
    assert sum((-1) ** y * (u - h) for y, u, h in degrees) == 0
    return degrees


def colour_singlets(table, length):
    weights = {v: n for (x, y, degree), n in table.items() if degree == length
               for v in signed_orbit(x, y)}
    return trivial_multiplicity(weights, "USp(4)")


def calculate():
    original = singlet_bound(AXIS_RECORDS)
    checks = check_controls()
    result = []
    shape = partner((7, 5, 4), 2)
    assert shape == (5, 5, 2, 1, 1)
    high, tail = shape[0], shape[2:]
    for length in (10, 12, 14):
        unprojected = original(shape, length)
        projected, _ = projection(high, tail, length, False)
        reduced, traces = projection(high, tail, length, True)
        assert 0 <= reduced <= projected
        assert page_trace(high, (1, 1, 1, 1), length, False, False) == original((5, 5, 1, 1, 1, 1), length)
        result.append(dict(source_length=length + 1, target_length=length,
                           exact_cochains=multiplicity(shape, length), prior_axis_bound=unprojected,
                           permutation_bound=projected, transverse_homology_bound=reduced,
                           cycle_class_traces=traces))
    return dict(group="USp(4)", source_flavour=[7, 5, 4], order=2, target_diagram=shape,
                four_slot_control=checks, results=result)
