"""Rational Levi and regular-element checks on the one-second-jet F4 source."""

import json
from itertools import combinations
from math import prod

from flint import fmpq_mat

from colour import Colour
from commuting_source import levi_strata, simple_coordinates
from levi_source_minors import branches


def evaluate():
    colour = Colour()
    roots = simple_coordinates(colour)
    simple = colour.data["raising"]
    cartan = fmpq_mat([colour.weights[i] for i in simple])
    profiles = levi_strata(colour)
    rows = []
    for profile in profiles["levi_subsystems"]:
        nodes = set(profile["simple_nodes"])
        complement = [i for i, root in enumerate(roots) if any(
            root[j] for j in range(4) if j not in nodes
        )]
        semisimple = cartan.solve(fmpq_mat([[int(i not in nodes)] for i in range(4)]))
        coefficients = dict(enumerate(semisimple.entries()))
        coefficients.update({simple[i]: 1 for i in nodes})
        action = fmpq_mat(52, 52)
        for i, coefficient in coefficients.items():
            for j, bracket in enumerate(colour.br[i]):
                for k, value in bracket.items():
                    action[k, j] += coefficient * value
        block = fmpq_mat([[action[i, j] for j in complement] for i in complement])
        determinant = prod(sum(roots[i][j] for j in range(4) if j not in nodes) for i in complement)
        assert (not complement or block.det() == determinant) and determinant
        assert all(not action[i, j] and not action[j, i]
                   for i in range(52) if i not in complement for j in complement)
        echelon, rank = action.rref()
        assert rank == 48
        pivots = [next(j for j in range(52) if echelon[i, j]) for i in range(rank)]
        centralizer = []
        for j in sorted(set(range(52)) - set(pivots)):
            column = {k: -echelon[i, j] for i, k in enumerate(pivots) if echelon[i, j]}
            centralizer.append({**column, j: 1})
        for first, second in combinations(centralizer, 2):
            bracket = [0] * 52
            for i, x in first.items():
                for j, y in second.items():
                    for k, v in colour.br[i][j].items():
                        bracket[k] += x * y * v
            assert not any(bracket)
        rank_one_factors = profile["positive_roots"] == len(nodes)
        minima = profile["required_nilpotent_slots"]
        if len(nodes) >= 2:
            assert min(minima) > profile["max_height"]
        rows.append(dict(nodes=sorted(nodes), centre=4-len(nodes),
                         scalar_shapes=profiles["scalar_coefficients"],
                         minimum_nilpotent_slots=minima,
                         maximum_positive_root_height=profile["max_height"],
                         complementary_dimension=len(complement),
                         complementary_determinant=str(determinant),
                         centralizer_dimension=52-rank, centralizer_abelian=True,
                         torus_or_product_of_rank_one_factors=rank_one_factors))
    assert sum(r["torus_or_product_of_rank_one_factors"] for r in rows) == 8
    return dict(ordinary_shape=[4, 4, 4, 2], scalar_coefficients=profiles["scalar_coefficients"],
                native_automorphism_checks=profiles["checked_chevalley_brackets"],
                strata=rows,
                semisimple_shapes_by_centre_dimension={str(z): branches(z) for z in (1, 2)})


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(evaluate(), indent=2))
