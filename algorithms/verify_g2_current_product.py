"""Evaluate the physical stress-current product and the G2 incoming line on one Cartan jet."""
import gzip
import json
from fractions import Fraction
from itertools import permutations
from pathlib import Path

from verify_g2_np_partner import HERE, canonical, fw, polynomial, run as verify_source
from fortuity.full.letters import Q_letter, is_ferm, letters_up_to


WEIGHTS = ((1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1), (0, 0))
FIRST = tuple(u for u, v in WEIGHTS)
SECOND = tuple(u - 2*v for u, v in WEIGHTS)
JET = {(0, 0, (0,)): FIRST, (0, 0, (1,)): SECOND,
       (1, 0, ()): FIRST, (0, 1, (2,)): SECOND}
ZERO = (0,) * 7


def trace(word):
    degree = sum(is_ferm(letter) for letter in word)
    if degree >= 2:
        return 0
    diagonal = [1] * 7
    for letter in word:
        diagonal = [x*y for x, y in zip(diagonal, JET.get(letter, ZERO))]
    return sum(diagonal)


def evaluate(expression):
    value = [Fraction(0), Fraction(0)]
    for multi, c in expression.items():
        degree = sum(is_ferm(l) for t in multi for l in t)
        if degree < 2:
            for t in multi:
                c *= trace(t)
            value[degree] += c
    return tuple(value)


def factors(flavour, direction):
    jet = (1, 0) if direction == 0 else (0, 1)
    l = canonical({(((0, 0, (flavour,)), jet + ((),)),): 1})
    bosons = {}
    for order in permutations(range(3)):
        if order[0] == flavour:
            sign = (-1)**sum(order[i] > order[j] for i in range(3) for j in range(i+1, 3))
            term = (((0, 0, (order[1],)), jet + ((order[2],),)),)
            bosons[term] = Fraction(sign, 2)
    complementary = tuple(i for i in range(3) if i != flavour)
    matter = {(((0, 0, complementary), jet + ((),)),): (-1)**flavour}
    g = canonical({**bosons, **matter})
    scalar_residual = canonical(fw.Q_expr(bosons))
    fermion_residual = canonical(fw.Q_expr(matter))
    assert len(scalar_residual) == 1
    assert scalar_residual == {w: -c for w, c in fermion_residual.items()}
    assert not canonical(fw.Q_expr(l)) and not canonical(fw.Q_expr(g))
    return l, g


def run():
    prior = verify_source()
    assert prior["source_cohomology_dimension"] == 1
    archive = json.loads(gzip.decompress((HERE.parent / "data/N4_G2_np_partner_replay_20261118.json.gz").read_bytes()))
    # On this jet the only odd value is a multiple of the same exterior generator.
    # Verify every matrix component of every classical letter variation.
    checked = 0
    for letter in letters_up_to(1, 1):
        for color in range(7):
            values = [0, 0]
            for c, x, y in Q_letter(letter):
                degree = is_ferm(x) + is_ferm(y)
                if degree < 2:
                    values[degree] += c * JET.get(x, ZERO)[color] * JET.get(y, ZERO)[color]
            assert values == [0, 0]
            checked += 1
    metric = [[sum(x*y for x, y in zip(a, b)) for a in (FIRST, SECOND)] for b in (FIRST, SECOND)]
    assert metric == [[4, 0], [0, 12]]
    result = {}
    for flavour in range(3):
        for direction, weight in ((0, Fraction(1, 6)), (1, -Fraction(1, 6))):
            left, _ = factors(flavour, direction)
            _, right = factors(flavour, 1-direction)
            for w, c in left.items():
                for v, h in right.items():
                    result[w+v] = result.get(w+v, 0) + weight*c*h
    result = canonical(result)
    assert result == polynomial(archive["representative"])
    left, _ = factors(0, 0)
    _, right = factors(0, 1)
    assert evaluate(left) == (0, 4) and evaluate(right) == (6, 0)
    assert evaluate(result) == (0, 4)
    # Evaluate (4.83)-(4.84) directly on diagonal fields. Every matter
    # fermion vanishes on this jet, so its component convention drops out.
    def trace_pair(a, b):
        return sum(i*j for i,j in zip(a,b))
    scalar = [FIRST, SECOND, ZERO]
    gaugino = [FIRST, ZERO]
    scalar_jet = [[ZERO, ZERO, ZERO],[ZERO, ZERO, SECOND]]
    physical = []
    for i, j, orientation in ((0,1,1),(1,0,-1)):
        for flavour in range(3):
            l_value = trace_pair(scalar[flavour],gaugino[i])
            n = (flavour+1)%3
            k = (flavour+2)%3
            g_value = Fraction(trace_pair(scalar[n],scalar_jet[j][k])
                               -trace_pair(scalar[k],scalar_jet[j][n]),2)
            physical.append(orientation*l_value*g_value/6)
    assert physical == [4,0,0,0,0,0] and sum(physical)==evaluate(result)[1]
    return dict(cartan_trace_form=metric, vanishing_Q0_letter_color_evaluations=checked,
                current_relative_sign_checks=6, reconstructed_trace_terms=len(result),
                current_L_exterior_coefficient=4, current_G_scalar_coefficient=6,
                saved_line_exterior_coefficient=4, perturbative_source_dimension=1,
                physical_singlet_cartan_terms=[int(c) for c in physical])


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
