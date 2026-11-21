"""Exact coordinate and supercharge identities for the two three-colour scalar charts."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import sympy as sp

from su3_square_zero_image import ALPH
from fortuity.full.letters import Q_letter, is_ferm, theta_sign


def superfield_rules():
    # Psi = A + theta1 B + theta2 C + theta1 theta2 D, with coefficients on the right.
    pieces = ((), (0,), (1,), (0, 1))
    rules = (
        ((0, 0, 1),),
        ((0, 1, 1), (1, 0, -1)),
        ((0, 2, 1), (2, 0, -1)),
        ((0, 3, 1), (3, 0, 1), (1, 2, 1), (2, 1, -1)),
    )
    count = 0
    for z, w, theta in ALPH:
        field = pieces.index(tuple(i for i in theta if i < 2))
        rest = tuple(i for i in theta if i >= 2)
        splits = Counter()
        for first, second, sign in rules[field]:
            for n in range(len(rest) + 1):
                for left in combinations(rest, n):
                    right = tuple(i for i in rest if i not in left)
                    for p in range(z + 1):
                        for m in range(w + 1):
                            u = (p, m, pieces[first] + left)
                            v = (z-p, w-m, pieces[second] + right)
                            if u in ALPH and v in ALPH:
                                phase = (-1)**(len(rest) + is_ferm(u)*len(right))
                                splits[u, v] += sign * phase * theta_sign(left, right)
        assert dict(splits) == {(u, v): c for c, u, v in Q_letter((z, w, theta))}
        count += len(splits)
    return count


def cyclic_pair():
    p, q, r, t = sp.symbols("p q r t")
    a, b, c, x, y, z = sp.symbols("a b c x y z")
    B = sp.Matrix([[0, a, b], [1, p, r], [0, q, -p]])
    C = sp.Matrix([[0, x, c], [0, y, t], [1, z, -y]])
    F = B*C - C*B
    slots = ((0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (1, 2))
    independent = sp.Matrix([F[i, j] for i, j in slots])
    coordinates = (x, y, z, b, a, c)
    determinant = sp.expand(independent.jacobian(coordinates).det())
    assert abs(determinant) == 1
    prescribed = sp.symbols("k0:6")
    inverse, = sp.solve(list(independent - sp.Matrix(prescribed)), coordinates, dict=True)
    assert all(sp.expand(v.subs(inverse) - s) == 0 for v, s in zip(independent, prescribed))
    zero = {v: val.subs(dict.fromkeys(prescribed, 0)) for v, val in inverse.items()}
    B, C = B.subs(zero), C.subs(zero)
    assert (B*C-C*B).applyfunc(sp.expand) == sp.zeros(3)
    traces = sp.Matrix([sp.trace(B**2), sp.trace(B**3), sp.trace(B*C), sp.trace(B**2*C)])
    jac = traces.jacobian((p, q, r, t))
    variable = sp.Symbol("X")
    disc = sp.discriminant(B.charpoly(variable).as_expr(), variable)
    assert sp.expand(jac.det(method="domain-ge") - 6*disc) == 0
    regular = {p: 0, q: 1, r: 1, t: 0}
    assert disc.subs(regular) == 108
    D = sp.Matrix([[0, x, y], [0, 0, 0], [0, 0, 0]])
    assert (sp.trace(D*B), sp.trace(D*C)) == (x, y)
    return dict(commutator_slots=[[i+1, j+1] for i, j in slots],
                transverse_jacobian=int(determinant), commuting_parameters=4,
                commuting_matrices=[str(B), str(C)], trace_jacobian_over_discriminant=6,
                regular_sample_discriminant=108, regular_sample_jacobian=int(jac.subs(regular).det()))


def subregular():
    p, q, u, r = sp.symbols("p q u r")
    a, b, c, d, e, f, g, h = sp.symbols("a b c d e f g h")
    B = sp.Matrix([[0, p, u], [1, q, 0], [0, r, -q]])
    C = sp.Matrix([[a, b, c], [d, e, f], [g, h, -a-e]])
    slots = ((0, 0), (1, 0), (2, 0), (1, 2))
    F = B*C - C*B
    eliminated = (b, e, h, c)
    independent = sp.Matrix([F[i, j] for i, j in slots])
    det = independent.jacobian(eliminated).det()
    assert abs(det) == 1
    substitution, = sp.solve(independent, eliminated, dict=True)
    A = C.subs(substitution)
    # Odd coefficients of the reduced A, in the ordered exterior basis a,d,f,g.
    ghosts = (a, d, f, g)
    generator = [A.diff(k) for k in ghosts]
    assert A == sum((k*M for k, M in zip(ghosts, generator)), sp.zeros(3)).applyfunc(sp.expand)
    even = (p, q, u, r)
    even_slots = ((0, 1), (1, 1), (0, 2), (2, 1))
    odd_slots = ((0, 0), (1, 0), (1, 2), (2, 0))
    QB = A*B-B*A
    variations = [sp.expand(QB[i, j]) for i, j in even_slots]
    assert (QB - sum((s*B.diff(k) for k, s in zip(even, variations)), sp.zeros(3))).applyfunc(sp.expand) == sp.zeros(3)
    odd_variations = [sp.S.Zero for _ in ghosts]
    for i, j in combinations(range(4), 2):
        # Compare the coefficient of each exterior product in Q(A) and A squared.
        square = generator[i]*generator[j] - generator[j]*generator[i]
        values = [square[m, n] for m, n in odd_slots]
        for k, value in enumerate(values):
            odd_variations[k] += ghosts[i]*ghosts[j]*value
        by_rule = sum((value*M for value, M in zip(values, generator)), sp.zeros(3))
        for k, qb in zip(even, variations):
            by_rule += qb.diff(ghosts[i])*generator[j].diff(k) - qb.diff(ghosts[j])*generator[i].diff(k)
        assert (square-by_rule).applyfunc(sp.expand) == sp.zeros(3)
    expected = (-2*q*f*g, f*g, 3*a*f, -3*a*g)
    assert all(sp.expand(x-y) == 0 for x, y in zip(odd_variations, expected))
    minors = sp.Matrix([[p-2*q**2, u, r], [3*a, f, g]])
    m12, m13, m23 = (minors[:, pair].det() for pair in combinations(range(3), 2))
    assert (F.subs(substitution)-sp.Matrix([[0, -2*q*m23, m12], [0, m23, 0], [0, -m13, -m23]])).applyfunc(sp.expand) == sp.zeros(3)
    return dict(transverse_jacobian=int(det), even_superfields=[str(v) for v in even],
                odd_superfields=[str(v) for v in ghosts], even_variations=list(map(str, variations)),
                odd_variations=list(map(str, expected)), checked_exterior_colour_coefficients=6*9,
                commuting_minors=list(map(str, (m12, m13, m23))))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires assertions")
    record = dict(colours=3, source_letters=len(ALPH), source_signed_splits=superfield_rules(),
                  cyclic_pair=cyclic_pair(), subregular=subregular())
    args.out.write_text(json.dumps(record, indent=2)+"\n")
