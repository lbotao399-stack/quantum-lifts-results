"""Project all four-matter-fermion colour contacts on the R-007 two-gaugino Cartan component."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from math import prod
import json

from cartan_contacts import (
    OFF_DIAGONAL,
    SCALAR_UNITS,
    cartan_scalar_expansion,
    radial_support,
)
from cubic_colours import source
from source_fluctuations import sign


def cyclic(word):
    return min(word[k:] + word[:k] for k in range(len(word)))


def canonical(traces):
    return tuple(sorted(cyclic(tuple(w)) for w in traces))


def contracts():
    # X,X,X,Y are the four derivative colours. U,V are the spin-paired
    # gaugino colours, and Z is the colour of the output matter fermion.
    labels = "XXXYUVZ"
    products = set()
    for permutation in permutations(range(7)):
        if any(j == i for i, j in enumerate(permutation)):
            continue
        visited = set()
        traces = []
        for start in range(7):
            if start in visited:
                continue
            seq = []
            i = start
            while i not in visited:
                visited.add(i)
                seq.append(labels[i])
                i = permutation[i]
            traces.append(seq)
        products.add(canonical(traces))
    return sorted(products)


def projected(pattern):
    orbit = Counter()
    for reverse, swap in product((False, True), repeat=2):
        tracing = []
        for tr in pattern:
            tr = tr[::-1] if reverse else tr
            tracing.append([{"U": "V", "V": "U"}.get(x, x) if swap else x for x in tr])
        orbit[canonical(tracing)] += (-1) ** (reverse + swap)
    return {p: coeff for p, coeff in orbit.items() if coeff}


def polarize(orbit):
    terms = Counter()
    for pattern, c in orbit.items():
        for perm in permutations(range(3)):
            axes = iter(perm)
            terms[
                canonical(
                    [
                        [next(axes) if x == "X" else "YUVZ".index(x) + 3 for x in tr]
                        for tr in pattern
                    ]
                )
            ] += c
    return {p: coeff for p, coeff in terms.items() if coeff}


# Matrix units are encoded by 3*row+column. The four diagonal matrices below
# are H0,H1 and three times their trace-duals. Division by nine restores both
# duals in each nonzero contraction.
DIAGONALS = {9: (1, 0, -1), 10: (0, 1, -1), 11: (2, -1, -1), 12: (-1, 2, -1)}


@cache
def colour_trace(word):
    off = [i for i, x in enumerate(word) if x < 9]
    if not off:
        return sum(prod(DIAGONALS[x][j] for x in word) for j in range(3))
    if len(off) == 1:
        return 0
    start = off[0]
    origin = word[start] // 3
    j = origin
    scalar = 1
    for x in word[start:] + word[:start]:
        if x >= 9:
            scalar *= DIAGONALS[x][j]
        else:
            row, col = divmod(x, 3)
            if j != row:
                return 0
            j = col
    return scalar if j == origin else 0


def configurations():
    # Only two off-diagonal fermions occur in the source. The unremoved
    # fermion is diagonal, and the two root directions must be opposite.
    for root_slots in combinations(range(4), 2):
        diagonal_slots = [i for i in range(4) if i not in root_slots]
        for p, q in permutations(range(3), 2):
            for first, second, z in product(range(2), repeat=3):
                directions = [None] * 4
                directions[root_slots[0]] = 3 * p + q
                directions[root_slots[1]] = 3 * q + p
                directions[diagonal_slots[0]] = 4 * first
                directions[diagonal_slots[1]] = 4 * second
                colours = [
                    11 + i // 4 if i in (0, 4) else 3 * (i % 3) + i // 3
                    for i in directions
                ]
                yield tuple(directions), tuple(colours + [9, 10, 9 + z]), z


def contraction(circuit, colours):
    return sum(
        weight * prod(colour_trace(tuple(colours[i] for i in tr)) for tr in traces)
        for traces, weight in circuit.items()
    )


def difference(vector, row, multiple):
    for key, value in row.items():
        updated = vector.get(key, 0) - multiple * value
        if updated:
            vector[key] = updated
        else:
            vector.pop(key, None)


def independent(vectors):
    basis = {}
    selected = []
    for i, v in enumerate(vectors):
        candidate = dict(v)
        for pivot, row in sorted(basis.items()):
            if pivot in candidate:
                difference(candidate, row, candidate[pivot])
        if candidate:
            pivot = min(candidate)
            scale = Fraction(1, candidate[pivot])
            basis[pivot] = {key: value * scale for key, value in candidate.items()}
            selected.append(i)
    return selected


def differentiate(source_by_mask, directions, z):
    result = Counter()
    for flavour in range(3):
        slots = [8 * i + j for i, j in zip((0, 1, 2, flavour), directions)]
        wanted = sum(1 << j for j in set(slots))
        if wanted.bit_count() != 4:
            continue
        marked = 1 << (8 * flavour + 4 * z)
        for mask, coefficients in source_by_mask.items():
            if mask & wanted != wanted:
                continue
            tail = mask
            parity = 1
            for slot in slots:
                bit = 1 << slot
                parity *= sign(bit, tail ^ bit)
                tail ^= bit
            if marked & tail:
                continue
            parity *= sign(marked, tail)
            for s, c in coefficients:
                result[s, marked | tail] += parity * c
    return {key: v for key, v in result.items() if v}


def radial_coordinates(diagonal):
    # Dot orders all first Cartan entries before all second entries. Translate
    # the middle public fermion with its minus sign. Keep U=H0 and V=H1.
    result = Counter()
    for (left, right), coefficient in diagonal.items():
        s = mask = 0
        for axis, group in enumerate((left, right)):
            for p, m, theta in group:
                if not theta:
                    if (p, m, axis) not in ((2, 0, 0), (0, 2, 1)):
                        coefficient = 0
                        break
                    bit = 1 << (24 + axis)
                elif theta.bit_count() == 1:
                    s += SCALAR_UNITS[2 * (theta.bit_length() - 1) + axis]
                    continue
                else:
                    assert theta in (3, 5, 6)
                    i = (7 ^ theta).bit_length() - 1
                    coefficient *= -1 if i == 1 else 1
                    bit = 1 << (8 * i + 4 * axis)
                assert not mask & bit
                coefficient *= sign(mask, bit)
                mask |= bit
        if coefficient:
            assert mask >> 24 == 3
            result[s, mask & ((1 << 24) - 1)] += coefficient
    return {k: v for k, v in result.items() if v}


def apply_row(row, derivatives):
    total = Counter()
    for i, c in row.items():
        for key, v in derivatives[i].items():
            total[key] += c * v / 9
    return {key: v for key, v in total.items() if v}


def separating_coordinates(images, target):
    # Retain the target value at every elimination step. Each returned
    # coordinate combination is zero on every image and equals one on T.
    pivots = {}
    for coordinate in sorted(set(target).union(*images)):
        column = {
            i: row[coordinate] for i, row in enumerate(images) if coordinate in row
        }
        functional = {coordinate: Fraction(1)}
        for pivot, (old_column, old_functional) in sorted(pivots.items()):
            if pivot in column:
                multiple = column[pivot]
                difference(column, old_column, multiple)
                difference(functional, old_functional, multiple)
        if column:
            pivot = min(column)
            scale = Fraction(1, column[pivot])
            pivots[pivot] = (
                {key: c * scale for key, c in column.items()},
                {key: c * scale for key, c in functional.items()},
            )
        else:
            value = sum(c * target.get(key, 0) for key, c in functional.items())
            if value:
                separator = {key: c / value for key, c in functional.items()}
                assert all(
                    sum(c * row.get(k, 0) for k, c in separator.items()) == 0
                    for row in images
                )
                return separator
    raise ValueError("The projected target lies in the image")


def show_coordinate(key, factor):
    scalar, mask = key
    return dict(
        scalar_powers=[scalar // unit % 6 for unit in SCALAR_UNITS],
        fermion_cartan_slots=[
            [i + 1, j + 1]
            for i, j in product(range(3), range(2))
            if mask & (1 << (8 * i + 4 * j))
        ],
        weight=str(factor),
    )


def determinant(entries):
    return sum(
        (-1) ** sum(i > j for a, i in enumerate(p) for j in p[a + 1 :])
        * prod(entries[i][j] for i, j in enumerate(p))
        for p in permutations(range(len(entries)))
    )


def projected_basis():
    terms = source()
    polynomial, support = cartan_scalar_expansion(terms)
    assert support == [0, 0, 3986, 0, 0, 0]
    source_by_mask = defaultdict(list)
    for (s, mask), c in polynomial.items():
        assert (mask & OFF_DIAGONAL).bit_count() == 2
        source_by_mask[mask].append((s, c))
    points = list(configurations())
    assert len(points) == 288
    patterns = contracts()
    circuits, labels = [], []
    for pattern in patterns:
        orbit = projected(pattern)
        if orbit and pattern == min(orbit):
            circuits.append(polarize(orbit))
            labels.append(pattern)
    tested = [
        {
            j: c
            for j, (_, colours, _) in enumerate(points)
            if (c := contraction(circuit, colours))
        }
        for circuit in circuits
    ]
    span = independent(tested)
    derivatives = [differentiate(source_by_mask, d, z) for d, _, z in points]
    images = [apply_row(tested[j], derivatives) for j in span]
    target = radial_coordinates(radial_support(terms)[1])
    image_span = independent(images)
    augmented_span = independent(images + [target])
    assert (
        len(patterns),
        len(circuits),
        len(span),
        len(image_span),
        len(augmented_span),
    ) == (337, 86, 3, 3, 4)
    assert len(target) == 80
    separator = separating_coordinates([images[i] for i in image_span], target)
    entries = [
        [row.get(key, 0) for row in [*(images[i] for i in image_span), target]]
        for key in sorted(separator)
    ]
    det = determinant(entries)
    assert det == 311115644928
    assert sum(c * target.get(key, 0) for key, c in separator.items()) == 1
    report = dict(
        trace_products=len(patterns),
        projected_trace_patterns=len(circuits),
        polarizations=len(points),
        colour_rank=len(span),
        source_image_rank=len(image_span),
        rank_with_radial_target=len(augmented_span),
        target_coordinates=len(target),
        independent_source_patterns=[
            ["".join(tr) for tr in labels[span[k]]] for k in image_span
        ],
        target_normalized_separator=[
            show_coordinate(key, c) for key, c in sorted(separator.items())
        ],
        distinguished_coefficients=[[str(v) for v in row] for row in entries],
        four_column_minor=str(det),
    )
    return report, [images[i] for i in image_span], target, separator


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(projected_basis()[0], indent=2))
