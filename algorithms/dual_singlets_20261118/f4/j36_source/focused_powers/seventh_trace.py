"""Identify the level-42 trace kernel with explicit products of level-36 cycles."""

from fractions import Fraction
import gzip
import json
from pathlib import Path
import sys

from sympy import Matrix

FIRST_PAGE = Path(__file__).resolve().parent.parent / "first_page"
sys.path.insert(0, str(FIRST_PAGE))
from trace_source_cycles import GLUE, canonical, coordinates


def basis(entries, columns):
    result = Matrix.zeros(len(entries["monomials"]), len(columns))
    for j, column in enumerate(columns):
        for i, value in column:
            result[i, j] = value
    return result


def run():
    sixth = json.loads(
        gzip.decompress((FIRST_PAGE / "trace_source_cycles.json.gz").read_bytes())
    )
    seventh = coordinates((5, 5, 2, 2, 2))
    target = {word: i for i, word in enumerate(seventh["monomials"])}
    product = Matrix.zeros(len(target), len(sixth["monomials"]))
    for j, raw in enumerate(sixth["monomials"]):
        sign, word = canonical((GLUE, *map(tuple, raw)))
        assert sign in (-1, 1)
        product[target[word], j] = sign
    assert product.T * product == Matrix.eye(162)

    ranks = {}
    for key in ("singlet_columns", "cycle_columns"):
        image = product * basis(sixth, sixth[key])
        rebuilt = basis(seventh, seventh[key])
        ranks[key] = [image.rank(), rebuilt.rank(), image.row_join(rebuilt).rank()]
    assert ranks == {"singlet_columns": [19, 19, 19], "cycle_columns": [14, 14, 14]}
    scores = json.loads((FIRST_PAGE / "trace_source_score.json").read_bytes())[
        "monomial_scores"
    ]
    source_scores = Matrix([[Fraction(x) for x in scores]])
    assert source_scores * basis(sixth, sixth["cycle_columns"]) == Matrix.zeros(1, 14)
    return {
        "dotted_flavour_weight": seventh["weight"],
        "reduced_trace_monomials": len(target),
        "monomials_by_number_of_glueball_traces": seventh["trace_power_distribution"],
        "classical_split_rank_before_quartic": seventh["tree_rank_before_quartic"],
        "classical_split_rank_after_quartic": seventh["tree_rank_after_quartic"],
        "ranks_product_direct_combined": ranks,
        "level_36_kernel_pairings": [0] * 14,
        "scope": "Reduced one-matter-fermion trace presentation with polarized quartic relations",
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
