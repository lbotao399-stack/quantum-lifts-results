"""Fixed nilpotent scalar fibres of the level-24 and level-30 one-fermion F4 sources."""
import argparse
from collections import Counter
from itertools import combinations, permutations
import json
from math import prod
from pathlib import Path


# Tables 8 and 3 of W. de Graaf, arXiv:math/0702193. Indices in the latter start at one.
POSITIVE_ROOTS = tuple(tuple(map(int, item)) for item in """
0001 1000 0010 0100 0011 1100 0110 0111
1110 0120 1111 0121 1120 1121 0122 1220
1221 1122 1231 1222 1232 1242 1342 2342
""".split())
ORBITS = {
    "1000": (24,), "0001": (21,), "0100": (17, 22), "2000": (16, 18),
    "0002": (11, 12), "0010": (14, 15, 16), "2001": (9, 15), "0101": (8, 14, 16),
    "1010": (10, 11, 15), "0200": (8, 9, 10, 18), "2200": (8, 2, 10),
    "1012": (10, 1, 9), "0202": (6, 5, 7, 13), "2202": (5, 4, 2, 10),
    "2222": (1, 3, 4, 2),
}
DOUBLE_SIMPLE = ((0, 2, -2, 0), (0, 0, 2, -2), (0, 0, 0, 2), (1, -1, -1, -1))


def add(left, right):
    return tuple(x+y for x, y in zip(left, right, strict=True))


def negative(weight):
    return tuple(-x for x in weight)


def minor(rows, cols):
    result = 0
    for order in permutations(cols):
        sign = (-1)**sum(i > j for k, i in enumerate(order) for j in order[k+1:])
        result += sign*prod(row[i] for row, i in zip(rows, order, strict=True))
    return result


def centralizer(diagram, indices):
    labels = tuple(map(int, diagram))
    first, second, third, fourth = labels
    h = (first+2*second+3*third+2*fourth, first+second+third, second+third, third)
    roots = {*POSITIVE_ROOTS, *map(negative, POSITIVE_ROOTS)}
    defining = tuple(POSITIVE_ROOTS[i-1] for i in indices)
    assert len(roots) == 48 and all(sum(x*y for x, y in zip(r, labels)) == 2 for r in defining)
    assert all(add(r, negative(s)) not in roots for r, s in permutations(defining, 2))
    assert any(minor(defining, cols) for cols in combinations(range(4), len(defining)))
    axes = tuple(combinations(range(4), len(defining)+1))
    euclidean = tuple(tuple(sum(c*simple[i] for c, simple in zip(r, DOUBLE_SIMPLE))
                            for i in range(4)) for r in defining)
    # h is in the span of the coroots. Opposite root vectors can therefore complete the triple.
    assert all(minor((*euclidean, h), cols) == 0 for cols in axes)

    def weight(root):
        return (sum(i*j for i, j in zip(root, labels)),
                *(minor((*defining, root), cols) for cols in axes))

    adjoint = Counter(weight(r) for r in roots)
    zero = weight((0, 0, 0, 0))
    adjoint[zero] += 4
    central = {}
    for (degree, *charge), multiplicity in adjoint.items():
        if degree >= 0:
            count = multiplicity-adjoint[(degree+2, *charge)]
            assert count >= 0
            if count:
                central[(degree, *charge)] = count
    assert sum((w[0]+1)*n for w, n in central.items()) == 52
    assert central.get(zero, 0) == 4-len(defining)
    # The torus perpendicular to the defining roots is a full Cartan of the reductive centralizer.
    for w, n in central.items():
        if w[0] == 0 and w != zero:
            assert n == central[negative(w)] == 1
    return central


def symmetric_characters(weights, zero, degree_limit, height_limit):
    series = [Counter() for _ in range(degree_limit+1)]
    series[0][zero] = 1
    for weight, multiplicity in weights.items():
        for _ in range(multiplicity):
            for k in range(1, degree_limit+1):
                for previous, count in series[k-1].items():
                    if previous[0]+weight[0] <= height_limit:
                        series[k][add(previous, weight)] += count
    return series


def weyl_denominator(central, zero):
    result = {zero: 1}
    for weight in central:
        if weight[0] == 0 and weight > zero:
            updated = Counter(result)
            for previous, c in result.items():
                updated[add(previous, weight)] -= c
            result = {w: c for w, c in updated.items() if c}
    return result


def convolution_at(left, right, wanted):
    return sum(c*right.get(add(wanted, negative(w)), 0) for w, c in left.items())


def fibre(central, width):
    maximum = max(w[0] for w in central)
    cap = maximum-2*(width-1)
    zero = (0,)*len(next(iter(central)))
    denominator = weyl_denominator(central, zero)
    symmetric = symmetric_characters(central, zero, width+1, cap)
    torus_count = invariant_count = 0
    for weight, copies in central.items():
        for shift, signed_multiplicity in denominator.items():
            wanted = add((weight[0]-2*(width-1), *weight[1:]), shift)
            schur = (convolution_at(symmetric[width], symmetric[width-1], wanted)
                     - convolution_at(symmetric[width+1], symmetric[width-2], wanted))
            assert schur >= 0
            if shift == zero:
                torus_count += copies*schur
            invariant_count += signed_multiplicity*copies*schur
    assert 0 <= invariant_count <= torus_count
    result = dict(width=width, maximum_scalar_height=cap, torus_invariants=torus_count,
                  reductive_invariants=invariant_count)
    bound = invariant_count
    if bound:
        degree_zero = {w: n for w, n in central.items() if w[0] == 0}
        assert degree_zero[zero] == 1 and len(degree_zero) == sum(degree_zero.values()) == 3
        minimum = min(w[0] for w in central if w[0] > 0)
        needed = (width-1)*minimum
        result.update(degree_zero_algebra="sl2", lowest_positive_scalar_height=minimum,
                      required_height_after_commutator=needed)
        # The scheme cut out by [Y_0,Z_0] consists of rank-at-most-one 3 x 2 matrices.
        if needed > cap:
            bound = 0
    result["cohomology_image_upper"] = bound
    return result


def calculate():
    rows = []
    for diagram, indices in ORBITS.items():
        central = centralizer(diagram, indices)
        height_profile = Counter()
        for weight, count in central.items():
            height_profile[weight[0]] += count
        rows.append(dict(diagram=diagram, representative_root_indices=indices,
                         orbit_dimension=52-sum(central.values()),
                         centralizer_by_height=dict(sorted(height_profile.items())),
                         sources={6*m: fibre(central, m) for m in (4, 5)}))
    assert len(rows) == 15
    assert all(not result["cohomology_image_upper"] for row in rows for result in row["sources"].values())
    return dict(source="de Graaf, arXiv:math/0702193, Tables 3 and 8",
                zero_scalar_image={24: 0, 30: 0}, nonzero_nilpotent_orbits=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Enable assertions to verify the root decompositions")
    args.out.write_text(json.dumps(calculate(), indent=2)+"\n")
