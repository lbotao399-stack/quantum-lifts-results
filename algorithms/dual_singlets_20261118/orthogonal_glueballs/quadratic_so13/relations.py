"""Multiply a determinant coefficient by its closed trace factors before projection."""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import factorial, prod
from auxiliary import SHAPES, quotient
from verify import encode, so_canon, PHI, GHOSTS, unit, scalar_projector
from orthogonal_glueball_relation_span import determinant
from rooted_trace import attach_rooted

ALPHA = PHI + GHOSTS


def norm_power(w):
    return prod(factorial(n) for n in Counter(w).values())


@lru_cache(None)
def strip_word(w, closed):
    # The specified product has all closed traces on the left.
    sign, sorted_factors = so_canon(tuple(tuple(ALPHA[a] for a in t) for t in closed))
    if not sign:
        return 0, ()
    remainder = list(w)
    for t in encode(sorted_factors):
        if t not in remainder:
            return 0, ()
        remainder.remove(t)
    merged = sorted_factors + tuple(tuple(ALPHA[a] for a in t) for t in remainder)
    sg, full = so_canon(merged)
    assert encode(full) == w
    return sign * sg, tuple(remainder)


@lru_cache(None)
def strip_sites(sites, closed):
    if not closed:
        return {sites: 1}
    w = closed[0]
    need = tuple(w.count(i) for i in (3, 4, 0, 1, 2))
    assert sum(need) == 2
    coefficient = 2 if w == (4, 3) else -2
    out = Counter()
    for j, old in enumerate(sites):
        if any(a < b for a, b in zip(old, need)):
            continue
        # Move the chosen odd variables from the front into their labelled planes.
        crossings = (need[0] + need[1]) * sum(v[0] + v[1] for v in sites[:j]) + need[
            1
        ] * (old[0] - need[0])
        new = tuple(a - b for a, b in zip(old, need))
        rest = sites[:j] + ((new,) if any(new) else ()) + sites[j + 1 :]
        for remaining, c in strip_sites(rest, closed[1:]).items():
            out[remaining] += (-1) ** crossings * coefficient * c
    return {s: c for s, c in out.items() if c}


def determinant_evaluation(k, point, native=False):
    shape = SHAPES[k][0]
    atoms, p, m, x = quotient(k)[1][point]
    inside = [
        (v, b, axes) for (closed, v), (b, axes) in zip(shape, atoms) if not closed
    ]
    external = [(v, axes) for (closed, v), (b, axes) in zip(shape, atoms) if closed]
    spec = dict(
        paths=[v for v, b, axes in inside],
        matrices=[unit(4, b) for v, b, axes in inside],
        scalar_directions=[[unit(3, a) for a in axes] for v, b, axes in inside],
        single=(p, m),
    )
    r = determinant(6, spec)
    r.project = scalar_projector(2, x)
    r.letters = 18 - sum(len(v) for v, axes in external)
    r.pf_degree = r.letters - sum(len(v) - 1 for v, b, axes in inside)
    assert r.pf_degree > 13
    if not native:
        attach_rooted(r)
    counts = Counter()
    for flavours in product(range(3), repeat=sum(len(a) for _, a in external)):
        it = iter(flavours)
        blocks = []
        degree = 0
        for v, axes in external:
            ax = iter(axes)
            word = []
            for label in v:
                if label == "s":
                    i = next(it)
                    degree += 1 << (3 * (3 * i + next(ax)))
                    word.append(i)
                else:
                    word.append(3 if label == "p" else 4)
            blocks.append(tuple(word))
        counts[tuple(blocks), degree] += 1
    if not external:
        return r.coefficient, r.cartan, r
    sign = (-1) ** (r.mask.bit_count() * (r.mask.bit_count() - 1) // 2)

    @lru_cache(None)
    def finish(items, degree, plane):
        weight = r.project.weights
        polynomial = r.product(items, plane)
        value = sum(
            c * weight.get(x + degree, 0)
            for (m, x, z), c in polynomial.items()
            if m == r.mask and z == r.used
        )
        return sign * value if plane else Fraction(sign * value, norm_power(items))

    def coefficient(w):
        total = 0
        for (blocks, degree), mult in counts.items():
            s, v = strip_word(w, blocks)
            if s:
                total += s * mult * finish(v, degree, False)
        return total

    def site(w):
        total = 0
        for (blocks, degree), mult in counts.items():
            for v, s in strip_sites(w, blocks).items():
                total += s * mult * finish(v, degree, True)
        return total

    return coefficient, site, r
