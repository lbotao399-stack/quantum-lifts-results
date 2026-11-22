"""Check the lifted-mode support of the charge-one cubic Dirac defect in SU(3)."""

import json
from itertools import permutations
from source_fluctuations import add, scale, wedge


def lifted(flavour, barred):
    return {1 << (6 * flavour + 4 + int(barred)): 1}


def antisymmetric_pair(a, b):
    return add(
        wedge(lifted(a, True), lifted(b, False)),
        scale(wedge(lifted(b, True), lifted(a, False)), -1),
    )


def derivative(poly, index):
    bit = 1 << index
    return {
        m ^ bit: (-1) ** ((m & (bit - 1)).bit_count()) * v
        for m, v in poly.items()
        if m & bit
    }


def charge(mask):
    return sum(
        bool(mask & (1 << (6 * a + 4))) - bool(mask & (1 << (6 * a + 5)))
        for a in range(4)
    )


def verify():
    # Common nonzero coefficient in Green--Kovacs (3.5)--(3.6) factors out.
    pieces = []
    for seq in permutations(range(4)):
        a, b, c, d = seq
        sign = (-1) ** sum(x > y for i, x in enumerate(seq) for y in seq[i + 1 :])
        pieces.append(
            scale(wedge(antisymmetric_pair(a, b), antisymmetric_pair(c, d)), sign)
        )
    action = add(*pieces)
    assert action and all(m.bit_count() == 4 and charge(m) == 0 for m in action)
    gradients = [derivative(action, i) for i in range(24)]
    for i, g in enumerate(gradients):
        if i % 6 < 4:
            assert not g
        else:
            assert g and all(m.bit_count() == 3 and abs(charge(m)) == 1 for m in g)

    # Residual colour matrices from the two lifted zero modes are off diagonal
    # between the instanton plane and the transverse colour. Test their full
    # matrix spaces against the embedded field-strength space.
    def unit(i, j):
        return [[int((a, b) == (i, j)) for b in range(3)] for a in range(3)]

    embedded = [unit(0, 1), unit(1, 0), [[1, 0, 0], [0, -1, 0], [0, 0, 0]]]
    off_diagonal = [
        unit(i, j) for i in range(3) for j in range(3) if (i == 2) != (j == 2)
    ]
    for f in embedded:
        for residual in off_diagonal:
            assert (
                sum(f[i][j] * residual[j][i] for i in range(3) for j in range(3)) == 0
            )
    return dict(
        quartic_terms=len(action),
        exact_derivatives_zero=16,
        lifted_derivatives_nonzero=8,
        cubic_terms_per_nonzero_derivative=[len(v) for v in gradients if v],
        pointwise_colour_pairings_zero=len(embedded) * len(off_diagonal),
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
