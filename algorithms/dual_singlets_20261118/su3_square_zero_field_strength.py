"""Reduce the rank-one scalar result by the exact components of the field-strength variation."""
import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp
from su3_square_zero_image import LETTER, SCALARS
from fortuity.full.letters import Q_letter


def vary_field_strength():
    result = []
    rules = {(LETTER[u], LETTER[v]): c for c,u,v in Q_letter((0,0,(0,1,2)))}
    for (scalar, fermion), sign in rules.items():
        if scalar in SCALARS:
            assert rules[fermion, scalar] == -sign
            result.append((scalar, fermion, sign))
    assert len(result) == 3 and len(rules) == 6
    return sorted(result)


def reduce_chart(terms, rules, pivot):
    scalars = sp.symbols("x1:4")
    coordinate = {s: x for (s,_,_), x in zip(rules, scalars, strict=True)}
    _, pivot_fermion, pivot_sign = rules[pivot]
    # For U=E12 the three entries 21,23,31 and the diagonal difference 22-11
    # vanish on its centralizer. Solve these four equations for the pivot flavour.
    def substitute(k, pos):
        if k != pivot_fermion or pos not in (3,4,5,6):
            return [(8*k+pos, 1)]
        terms = [(8*k, 1)] if pos == 4 else []
        for i, (_, psi, sign) in enumerate(rules):
            if i != pivot:
                coefficient = -sign * scalars[i] / (pivot_sign * scalars[pivot])
                terms.append((8*psi+pos, coefficient))
                if pos == 4:
                    terms.append((8*psi, -coefficient))
        return terms

    residual = defaultdict(int)
    for term in terms:
        coeff = Fraction(term["coefficient"])
        powers = []
        for k,pos in term["even"]:
            assert k in SCALARS and pos == 1
            powers.append(coordinate[k])
        row = {0: sp.Rational(coeff.numerator, coeff.denominator) * sp.prod(powers)}
        for k,pos in term["odd"]:
            following = defaultdict(int)
            for var, factor in substitute(k, pos):
                for mask, value in row.items():
                    if not mask & (1 << var):
                        sign = (-1)**((mask >> (var+1)).bit_count())
                        following[mask | (1 << var)] += sign * factor * value
            row = following
        for mask, value in row.items():
            residual[mask] += value
    nonzero = {mask:sp.cancel(v) for mask,v in residual.items() if sp.cancel(v) != 0}
    return dict(nonzero_monomials=len(nonzero),
                exterior_coefficients=[dict(odd=[[i//8,i%8] for i in range(mask.bit_length()) if mask & (1 << i)],
                                            coefficient=str(v)) for mask,v in sorted(nonzero.items())])


def primitive(terms, rules):
    """Replace one scalar times a fermion entry by f_13 or f_32 and solve over Q."""
    target = {(tuple(k for k,_ in row["even"]), tuple(map(tuple, row["odd"]))):
              sp.Rational(row["coefficient"]) for row in terms}
    candidates = set()
    for even, odd in target:
        for scalar, fermion, _ in rules:
            if scalar not in even:
                continue
            for pos in (5, 6):
                if (fermion,pos) in odd:
                    rest = list(even)
                    rest.remove(scalar)
                    candidates.add((tuple(rest), pos, tuple(v for v in odd if v != (fermion,pos))))
    candidates = sorted(candidates)
    columns = []
    for even, pos, odd in candidates:
        col = {}
        for scalar, fermion, sign in rules:
            inserted = (fermion,pos)
            if inserted in odd:
                continue
            exterior_sign = (-1)**sum(v < inserted for v in odd)
            key = (tuple(sorted(even + (scalar,))), tuple(sorted((inserted,) + odd)))
            col[key] = sign * (1 if pos == 5 else -1) * exterior_sign
        columns.append(col)
    rows = sorted(set(target).union(*(set(col) for col in columns)))
    matrix = sp.Matrix([[col.get(row,0) for col in columns] for row in rows])
    rhs = sp.Matrix([-target.get(row,0) for row in rows])
    solution, free = matrix.gauss_jordan_solve(rhs)
    solution = solution.subs({v:0 for v in free})
    assert matrix * solution == rhs
    result = []
    for (even, pos, odd), value in zip(candidates, solution, strict=True):
        if value:
            result.append(dict(scalars=even, field_strength_entry=[1,3] if pos==5 else [3,2],
                               odd=odd, coefficient=str(value)))
    return dict(equations=matrix.rows, candidate_columns=matrix.cols, rank=matrix.rank(),
                correction=result, residual_terms=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    data = json.loads(args.input.read_bytes())
    assert data["level"] == 34 and data["direct_image"] == "Q1"
    terms = data["gauginos_in_scalar_centralizer"]["one_independent_scalar"]["polynomial"]
    rules = vary_field_strength()
    result = dict(level=34, scalar_matrix="E12", input_terms=len(terms),
                  scalar_fermion_commutators=rules, charts={str(i+1): reduce_chart(terms,rules,i) for i in range(3)},
                  polynomial_correction=primitive(terms, rules))
    args.out.write_text(json.dumps(result, indent=2)+"\n")
