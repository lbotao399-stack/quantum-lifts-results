"""Test all flavour contractions for the mixed colour tensor against the four Ward boundaries."""

from collections import Counter
from fractions import Fraction
from itertools import permutations, product
import json
from pathlib import Path

from cartan_contacts import SCALAR_UNITS, cartan_scalar_expansion
from contact_ward import (
    LAYOUTS,
    boundary,
    constraint,
    field_strength_variation,
    primitive_relation,
)
from cubic_colours import source
from four_fermion_contacts import (
    contracts,
    determinant,
    independent,
    polarize,
    projected,
    projected_basis,
)
from mixed_contact import U, V, root_pair
from source_fluctuations import sign

HERE = Path(__file__).resolve().parent


def epsilon(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    return (-1) ** (int(i > j) + int(i > k) + int(j > k))


def images(terms):
    # i is the root fermion, j differentiates along U or V, k along
    # the output scalar. The scalar derivative has flavour h, the output s.
    output = [Counter() for _ in range(3)]
    for i, h in product(range(3), repeat=2):
        for p, q in permutations(range(3), 2):
            residual = root_pair(terms, i, h, p, q)
            if not residual:
                continue
            for j, k, s in product(range(3), repeat=3):
                weights = (
                    int(i == h) * epsilon(j, k, s),
                    int(j == h) * epsilon(k, i, s),
                    int(k == h) * epsilon(i, j, s),
                )
                assert sum(weights) == int(s == h) * epsilon(i, j, k)
                if not any(weights):
                    continue
                for a, b in product(range(2), repeat=2):
                    first, second = 1 << (8 * j + 4 * a), 1 << (8 * k + 4 * b)
                    if first == second:
                        continue
                    colour = (U[q] - U[p]) * V[a] - (V[q] - V[p]) * U[a]
                    for (scalar, mask), value in residual.items():
                        if mask & first and mask & second:
                            tail = mask ^ first
                            parity = sign(first, tail) * sign(second, tail ^ second)
                            key = scalar + SCALAR_UNITS[2 * s + b], tail ^ second
                            for n, w in enumerate(weights):
                                if w:
                                    output[n][key] += value * colour * parity * w
    return [{key: c for key, c in row.items() if c} for row in output]


def verify():
    archive = json.loads((HERE / "contact_ward.json").read_text())
    _, pure, target, functional = projected_basis()
    probes = sorted(functional)
    old = [{int(i): Fraction(c) for i, c in r.items()} for r in archive["rows"]]
    variation = field_strength_variation()
    terms = [(Fraction(c), w) for c, w in source()]
    mixed = images(terms)
    controls = [[v.get(k, 0) for v in mixed] for k in probes]
    for layout in LAYOUTS:
        vecs = images(boundary(layout, variation))
        controls.append(
            [sum(c * v.get(k, 0) for k, c in functional.items()) for v in vecs]
        )
    assert all(sum(v) == r.get(86, 0) for v, r in zip(controls, old))
    choices = independent(
        [{i: r.get(k, 0) for i, r in enumerate(pure) if k in r} for k in probes]
    )
    sub = [[row.get(probes[i], 0) for row in pure] for i in choices]
    base_det = determinant(sub)
    assert len(choices) == 3 and base_det
    pure_weights = []
    for j in range(86):
        coords = []
        for col in range(3):
            matrix = [
                row[:col] + [old[i].get(j, 0)] + row[col + 1 :]
                for row, i in zip(sub, choices)
            ]
            coords.append(Fraction(determinant(matrix), base_det))
        pure_weights.append(coords)
    equations, values = [], []
    for key in sorted(set(target).union(*pure, *mixed)):
        row = [
            sum(c * v.get(key, 0) for c, v in zip(coords, pure))
            for coords in pure_weights
        ]
        row += [v.get(key, 0) for v in mixed]
        equations.append({j: c for j, c in enumerate(row) if c})
        values.append(target.get(key, 0))
    source_rank = len(independent(equations))
    for old_row, extra in zip(old[4:], controls[4:]):
        equations.append(
            {i: c for i, c in old_row.items() if i < 86}
            | {86 + i: c for i, c in enumerate(extra) if c}
        )
        values.append(0)
    relation = archive["primitive_left_relation"]
    tests = [sum(c * v[j] for c, v in zip(relation, controls)) for j in range(3)]
    assert sum(tests) == 0
    selected = independent(equations)
    assert source_rank == 4 and len(selected) == 8
    circuits = [
        polarize(projected(p))
        for p in contracts()
        if projected(p) and min(projected(p)) == p
    ]
    additional = ("0f0ad0bce", "a0fb/c0d0e")
    assert len(circuits) == 86
    for layout in additional:
        terms = boundary(layout, variation)
        polynomial, _ = cartan_scalar_expansion(terms)
        pure_row = constraint(polynomial, {}, functional, circuits)
        mixed_columns = [
            sum(c * v.get(k, 0) for k, c in functional.items()) for v in images(terms)
        ]
        equations.append(
            pure_row | {86 + i: c for i, c in enumerate(mixed_columns) if c}
        )
        values.append(0)
    augmented_rows = [r | {89: c} if c else r for r, c in zip(equations, values)]
    selected, augmented = independent(equations), independent(augmented_rows)
    assert len(selected) == 9 and len(augmented) == 10
    small = [equations[i] for i in augmented]
    rhs = [values[i] for i in augmented]
    relation_full = primitive_relation(small)
    assert all(
        sum(c * r.get(j, 0) for c, r in zip(relation_full, small)) == 0
        for j in range(89)
    )
    discrepancy = sum(c * v for c, v in zip(relation_full, rhs))
    columns = independent(
        [{i: r[j] for i, r in enumerate(small) if j in r} for j in range(89)]
    )
    assert len(columns) == 9
    minor = [[r.get(j, 0) for j in columns] + [v] for r, v in zip(small, rhs)]
    # Fraction-free elimination from SymPy avoids a factorial determinant expansion.
    from sympy import Matrix

    det = Matrix(minor).det()
    assert discrepancy and det
    return dict(
        flavour_colours=3,
        joint_columns=89,
        source_coordinates=len(equations) - 6,
        source_rank=source_rank,
        joint_rank=len(selected),
        augmented_rank=len(augmented),
        additional_flavour_values_on_left_relation=list(map(str, tests)),
        new_values_on_eight_ward_rows=[[str(c) for c in row] for row in controls],
        independent_equations=selected,
        augmented_independent_equations=augmented,
        boundary_layouts=[*LAYOUTS, *additional],
        pivot_columns=columns,
        augmented_minor=[[str(v) for v in row] for row in minor],
        augmented_determinant=str(det),
        primitive_left_relation=relation_full,
        left_relation_on_target=str(discrepancy),
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
