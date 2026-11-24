"""Check the entire degree-three invariant bilinear complex without a subspace projector.

All 60 columns pair a first Taylor jet with a second Taylor jet using the
native invariant form.  The classical outputs use the nonzero alternating
native three-colour tensor.  The first quantum action is computed from the
native two-input Taylor rule and the complete adjoint Casimir contraction.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import sys

from flint import fmpz_mat


def letters(degree):
    return [(a, b, flavours) for a in range(degree + 1) for b in range(degree - a + 1)
            for flavours in combinations(range(3), degree - a - b)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    commit = "c778cd7c3ceb72ad8082d791b2aee5f0df082164"
    root = args.inputs / commit
    page = root / "algorithms/dual_singlets_20261118/f4/j36_source/first_page"
    runtime = root / "algorithms/su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work"
    sys.path[:0] = [str(page), str(runtime / "so8-bmn/vendor"), str(runtime / "so8-protected")]
    from colour import Colour
    from fortuity.full.letters import Q_letter, is_ferm
    from q1d_exact import pair_terms
    first, second = letters(1), letters(2)
    assert len(first) == 5 and len(second) == 12
    columns = [(marked, scalar) for marked in second for scalar in first]
    rows = [indices for indices in combinations_with_replacement(range(5), 3)
            if all(indices.count(i) <= 1 or is_ferm(first[i]) for i in set(indices))]
    assert len(columns) == 60 and len(rows) == 20
    row_index = {row: i for i, row in enumerate(rows)}

    def canonical(triple):
        indices = [first.index(letter) for letter in triple]
        if any(indices.count(i) > 1 and not is_ferm(first[i]) for i in indices):
            return None, 0
        sign = 1
        for i, a in enumerate(indices):
            for b in indices[i + 1:]:
                if a > b:
                    sign *= -(-1) ** (int(is_ferm(first[a])) * int(is_ferm(first[b])))
        return tuple(sorted(indices)), sign

    tree = [[0] * 60 for row in rows]
    for column, (marked, scalar) in enumerate(columns):
        assert not Q_letter(scalar)
        for coefficient, left, right in Q_letter(marked):
            target, sign = canonical((left, right, scalar))
            if sign:
                tree[row_index[target]][column] += int(coefficient) * sign

    colour = Colour()
    casimir = [[Fraction(0) for j in range(52)] for i in range(52)]
    for (a, b), inverse in colour.inverse.items():
        inverse = Fraction(str(inverse))
        for source in range(52):
            for middle, value in colour.br[b][source].items():
                for target, other in colour.br[a][middle].items():
                    casimir[target][source] += inverse * value * other
    assert casimir == [[Fraction(18 * (i == j)) for j in range(52)] for i in range(52)]
    ghosts = ((1, 0, ()), (0, 1, ()))
    quantum = []
    for marked, scalar in columns:
        terms = Counter()
        for coefficient, left, right in pair_terms(marked, scalar):
            terms[left, right] += Fraction(coefficient)
        assert all(pair in (ghosts, ghosts[::-1]) for pair, coefficient in terms.items() if coefficient)
        value = 18 * (terms[ghosts[::-1]] - terms[ghosts])
        assert value.denominator == 1
        quantum.append(int(value))
    scalar_triplet = tuple(first.index((0, 0, (i,))) for i in range(3))
    comparison = tree[row_index[scalar_triplet]]
    assert quantum == [18 * value for value in comparison]
    matrix = fmpz_mat(tree)
    assert matrix.rank() == 20
    kernel, nullity = matrix.nullspace()
    assert nullity == 40 and (matrix * kernel).is_zero()
    assert (fmpz_mat([quantum]) * kernel).is_zero()
    live = [i for i, value in enumerate(quantum) if value]
    assert len(live) == 3
    result = {
        "schema": "f4-complete-closed-bilinear-v1", "field": "Q", "public_commit": commit,
        "first_jets": first, "second_jets": second,
        "ordered_bilinear_columns": columns, "alternating_three_colour_rows": rows,
        "twice_classical_matrix": tree, "first_quantum_row": quantum,
        "nonzero_quantum_columns": live,
        "adjoint_casimir_on_every_colour": 18,
        "verified_casimir_entries": 52 * 52,
        "classical_rank": 20, "complete_closed_bilinear_dimension": 40,
        "first_quantum_rank_on_complete_classical_kernel": 0,
        "exact_row_identity": "first_quantum_row = 18 * twice_classical_matrix[three_distinct_scalar_flavours]",
        "all_higher_bare_operations_vanish": "For order at least two the bare operation requires at least three input letters; every bilinear has two.",
        "closed_first_jet_spectators": "Any order-r positive bare operation needs 3r flavour slots on r+1 letters; r+1 first jets supply at most r+1 flavour slots, which is less than 3r for every r at least one.",
        "scope": "Every classical closed invariant bilinear with three Taylor slots is strictly closed under the bare perturbative differential. Matched quantum multiplication with strictly closed first-jet invariant polynomials gives a formally closed completion of the exceptional product. No Taylor or colour subspace projector is used in this calculation.",
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in
                     ("first_jets", "second_jets", "ordered_bilinear_columns",
                      "alternating_three_colour_rows", "twice_classical_matrix", "first_quantum_row")}), flush=True)


if __name__ == "__main__":
    main()
