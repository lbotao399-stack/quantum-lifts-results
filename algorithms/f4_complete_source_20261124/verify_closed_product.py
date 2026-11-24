"""Independently solve every rational raising and marked-factor closure equation."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
import json
from math import gcd, lcm
from pathlib import Path
import sys


def parity(indices):
    return (-1) ** sum(a > b for j, a in enumerate(indices) for b in indices[j + 1:])


def equations():
    triples = list(combinations_with_replacement(range(4), 3))
    bilinears = list(combinations_with_replacement(range(4), 2))
    pairs = list(combinations(range(4), 2))
    columns = []
    for omitted, gram, mark, last in product(triples, bilinears, pairs, range(4)):
        weights = [3 - omitted.count(i) + gram.count(i) + mark.count(i) + (last == i) for i in range(4)]
        if weights == [4, 4, 4, 2]:
            columns.append((omitted, gram, mark, last))
    assert len(columns) == 96
    rows = defaultdict(lambda: [0] * 96)
    for position, (omitted, gram, mark, last) in enumerate(columns):
        for i in range(3):
            j = i + 1
            # The omitted index labels the dual of a three-fold alternating tensor.
            for at, index in enumerate(omitted):
                if index == i:
                    changed = list(omitted); changed[at] = j
                    rows[("raising", i, tuple(sorted(changed)), gram, mark, last)][position] -= 1
            for at, index in enumerate(gram):
                if index == j:
                    changed = list(gram); changed[at] = i
                    rows[("raising", i, omitted, tuple(sorted(changed)), mark, last)][position] += 1
            for at, index in enumerate(mark):
                if index == j:
                    changed = list(mark); changed[at] = i
                    if len(set(changed)) == 2:
                        rows[("raising", i, omitted, gram, tuple(sorted(changed)), last)][position] += parity(changed)
            if last == j:
                rows[("raising", i, omitted, gram, mark, i)][position] += 1
        labels = (*mark, last)
        if len(set(labels)) == 3:
            rows[("classical_split", omitted, gram, tuple(sorted(labels)))][position] += parity(labels)
    return columns, dict(rows)


def rational_kernel(integer_rows):
    matrix = [[Fraction(v) for v in row] for row in integer_rows]
    pivots = []
    for column in range(96):
        number = len(pivots)
        found = next((i for i in range(number, len(matrix)) if matrix[i][column]), None)
        if found is None:
            continue
        matrix[number], matrix[found] = matrix[found], matrix[number]
        pivot = matrix[number][column]
        matrix[number] = [v / pivot for v in matrix[number]]
        for i in range(len(matrix)):
            if i != number and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [a - factor * b for a, b in zip(matrix[i], matrix[number])]
        pivots.append(column)
    free = [i for i in range(96) if i not in pivots]
    assert len(pivots) == 95 and len(free) == 1
    vector = [Fraction(0)] * 96
    vector[free[0]] = Fraction(1)
    for i, pivot in enumerate(pivots):
        vector[pivot] = -matrix[i][free[0]]
    denominator = lcm(*(v.denominator for v in vector))
    return [int(v * denominator) for v in vector], pivots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    columns, rows = equations()
    coefficient, pivots = rational_kernel(list(rows.values()))
    residual = [sum(a * b for a, b in zip(row, coefficient)) for row in rows.values()]
    assert not any(residual) and sum(bool(a) for a in coefficient) == 87
    assert gcd(*coefficient) == 1
    sys.path.insert(0, str(args.inputs))
    from albert_product_cycles import highest
    reference_columns, reference = highest()
    assert columns == reference_columns and coefficient == list(map(int, reference[0]))
    spectators = sorted({(omitted, gram) for a, (omitted, gram, mark, last)
                         in zip(coefficient, columns) if a})
    # These equations hold for each fixed pair of spectator polynomials.
    # Their zero residual proves closure of the remaining marked bilinear.
    split_count = sum(key[0] == "classical_split" for key in rows)
    terms = [{"coefficient": a, "omitted_flavours": omitted, "scalar_bilinear": gram,
              "marked_pair": mark, "marked_bilinear_last": last}
             for a, (omitted, gram, mark, last) in zip(coefficient, columns) if a]
    result = {"field": "Q", "formal_columns": 96, "equations": len(rows),
              "rational_rank": len(pivots), "rational_kernel_dimension": 1,
              "nonzero_coefficients": 87, "coefficient_gcd": 1,
              "raising_equations": len(rows) - split_count,
              "classical_split_equations": split_count,
              "separately_closed_spectator_pairs": len(spectators),
              "maximum_absolute_constraint_residual": max(map(abs, residual)),
              "matches_public_rational_product": True, "terms": terms,
              "constraint_matrix_sha256": sha256(json.dumps(list(rows.values()), separators=(",", ":")).encode()).hexdigest(),
              "scope": "Exact rational closure and highest-weight equations. Evaluation minors separately certify native finite-colour independence."}
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "terms"}), flush=True)


if __name__ == "__main__":
    main()
