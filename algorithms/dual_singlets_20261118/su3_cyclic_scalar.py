"""Verify the companion coordinates and the classical superfield signs at three colours."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import sympy as sp
from su3_square_zero_image import ALPH
from fortuity.full.letters import Q_letter, is_ferm, theta_sign


def field_splits(letter):
    z, w, theta = letter
    signed = 0 in theta
    remaining = tuple(i for i in theta if i != 0)
    terms = Counter()
    # In Psi = A + theta^1 B the equations are Q0 A = A^2, Q0 B = A B - B A.
    products = [(False, True, 1), (True, False, -1)] if signed else [(False, False, 1)]
    for first_b, second_b, sign in products:
        for n in range(len(remaining) + 1):
            for first in combinations(remaining, n):
                second = tuple(i for i in remaining if i not in first)
                for p in range(z + 1):
                    for m in range(w + 1):
                        u = (p, m, ((0,) if first_b else ()) + first)
                        v = (z-p, w-m, ((0,) if second_b else ()) + second)
                        if u in ALPH and v in ALPH:
                            phase = (-1)**(len(remaining) + is_ferm(u) * len(second))
                            terms[u, v] += sign * phase * theta_sign(first, second)
    return terms


def run():
    split_count = 0
    for letter in ALPH:
        terms = field_splits(letter)
        original = {(u,v): c for c,u,v in Q_letter(letter)}
        assert dict(terms) == original
        split_count += len(terms)

    entries = sp.symbols("b0:8")
    matrix = sp.Matrix(3, 3, [*entries, -entries[0]-entries[4]])
    first = sp.eye(3)[:, 0]
    second = matrix * first
    third = matrix * second
    frame = first.row_join(second).row_join(third)
    quadratic, cubic = sp.trace(matrix**2)/2, matrix.det()
    companion = sp.Matrix([[0,0,cubic], [1,0,quadratic], [0,1,0]])
    assert (matrix * frame - frame * companion).applyfunc(sp.expand) == sp.zeros(3)
    coordinates = second.col_join(third).col_join(sp.Matrix([quadratic,cubic]))
    jacobian = coordinates.jacobian(entries).to_DM()
    determinant = jacobian.domain.to_sympy(jacobian.det())
    assert sp.expand(determinant-frame.det()**2) == 0
    principal = sp.Matrix([[0,0,0],[1,0,0],[0,1,0]])
    nil_frame = first.row_join(principal*first).row_join(principal**2*first)
    assert principal**3 == sp.zeros(3) and principal**2 != sp.zeros(3) and nil_frame.det() == 1
    return dict(colours=3, letters_checked=len(ALPH), signed_splits_checked=split_count,
                generic_companion_residual_entries=0, jacobian_over_cyclic_determinant_squared="1",
                principal_nilpotent_cyclic_determinant=int(nil_frame.det()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    args.out.write_text(json.dumps(run(), indent=2)+"\n")
