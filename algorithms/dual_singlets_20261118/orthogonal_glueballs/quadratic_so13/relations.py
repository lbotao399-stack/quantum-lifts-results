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


def merge_double_paths(relation, indices):
    """Extract a factorial times one commuting power for identical even paths."""
    count = len(indices)
    if count < 2:
        return
    mask = sum(1 << i for i in indices)
    first = 1 << indices[0]
    assert relation.used & mask == mask
    project = relation.project
    stride = 1 << (max(project.divisors).bit_length() + 1)
    for block, choices in relation.catalog.items():
        representative = {(b, m, x, c) for b, m, x, t, c in choices if t == first}
        for i in indices[1:]:
            assert representative == {
                (b, m, x, c) for b, m, x, t, c in choices if t == (1 << i)
            }
        assert all(m == 0 for b, m, x, c in representative)
        relation.catalog[block] = [v for v in choices if not (v[3] & mask)]
        relation.catalog[block].extend(
            (b, m, x + stride, 0, c) for b, m, x, c in representative
        )
    relation.used ^= mask
    project.weights = {
        x + count * stride: factorial(count) * c for x, c in project.weights.items()
    }
    project.divisors = {
        x + j * stride for x in project.divisors for j in range(count + 1)
    }


def determinant_evaluation(
    k,
    point,
    native=False,
    *,
    shapes=SHAPES,
    length=18,
    colours=13,
    compress_pairs=False,
    extra_axes=2,
):
    shape = shapes[k][0]
    atoms, p, m, x = quotient(k, shapes, extra_axes)[1][point]
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
    r.letters = length - sum(len(v) for v, axes in external)
    r.pf_degree = r.letters - sum(len(v) - 1 for v, b, axes in inside)
    assert r.pf_degree > colours
    if compress_pairs:
        merge_double_paths(
            r, [i for i, (v, b, _) in enumerate(inside) if v == "mp" and b == 3]
        )
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
