"""Compare ten-letter F4 trace relations with native Albert ghost-bracket boundaries."""
import argparse
from functools import lru_cache
from itertools import permutations
import json
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat, nmod_mat

from jordan_bracket_singlets import LETTERS, degree, singlets
from native_trace_relations import full_reduction, obstruction_basis, terminal_dimension
from scalar_orthogonal_walls import native, reduce_native
from trace_source_cycles import quartic_pivots
from f4_colour import basis as albert_basis, jordan_colours
from f4_subgroup_equation import integral_rows
from fortuity.hwcoh.evalg import table


def maximum(a):
    return int(max(abs(int(a.max())), abs(int(a.min()))))


def checked(a, b, factor):
    # Include the entire exterior sum before any reduction or contraction.
    if maximum(a) * maximum(b) * factor >= 2 ** 62:
        return a.astype(object), b.astype(object)
    return a, b


def contract(frame, graphs, start, end):
    left = jordan_colours()
    gram = np.einsum("aij,bji->ab", left, left)
    metric = [int(x) for x in gram.diagonal()]
    scale = lcm(*metric)
    inverse = np.array([scale // x for x in metric], dtype=np.int64)
    jordan_duals = left.transpose(0, 2, 1).reshape(26, -1).T
    k = frame.k

    def multiply(a, da, b, db):
        i, j, signs, starts = table(k, da, db)
        x, y = checked(a, b, 27 * comb(da + db, da))
        terms = (x[:, i] @ y[:, j]) * signs[None, :, None, None]
        return np.add.reduceat(terms, starts, axis=1)

    @lru_cache(None)
    def lie(node):
        if node[0] == 0:
            return frame.L[LETTERS[node[1]]][1][start:end]
        _, x, y = node
        a, b = lie(x), lie(y)
        da, db = degree(x), degree(y)
        ab, ba = multiply(a, da, b, db), multiply(b, db, a, da)
        return ab - (-1) ** (da * db) * ba

    @lru_cache(96)
    def word(slots):
        if len(slots) == 1:
            return lie(slots[0])
        return multiply(word(slots[:-1]), sum(degree(s) for s in slots[:-1]),
                        lie(slots[-1]), degree(slots[-1]))

    @lru_cache(None)
    def covariant(arm):
        matrices = word(arm[1:])
        a, b = checked(matrices, jordan_duals, 27 ** 2)
        return a.reshape(*a.shape[:2], -1) @ b

    def pair(a, b):
        da, db = degree(a), degree(b)
        i, j, signs, starts = table(k, da, db)
        v, w = checked(covariant(a), covariant(b),
                       26 * int(inverse.max()) * comb(da + db, da))
        terms = (v[:, i] * w[:, j] * inverse).sum(axis=-1) * signs
        return np.add.reduceat(terms, starts, axis=1)

    values = []
    for a, b, tail in graphs:
        coeff = pair(a, b)
        if tail[1:]:
            mat = word(tail[1:])
            if maximum(mat) * 27 >= 2 ** 62:
                mat = mat.astype(object)
            scalar = mat.trace(axis1=-2, axis2=-1)
            deg = degree(tail)
            dx = degree(a) + degree(b)
            x, y = checked(coeff, scalar, comb(deg + dx, dx))
            coeff, _ = frame.scal_mul(x, dx, y, deg)
        values.append(coeff.reshape(-1))
    # Verify the graded quartic symmetry before contracting the two arms.
    for a in sorted({arm for graph in graphs for arm in graph[:2]})[::43]:
        for perm in permutations(range(4)):
            slots = a[1:]
            odd = [i for i in perm if degree(slots[i]) % 2]
            phase = (-1) ** sum(x > y for i, x in enumerate(odd) for y in odd[i+1:])
            alternative = (2, *(slots[i] for i in perm))
            assert np.array_equal(covariant(a), phase * covariant(alternative))
    return np.array(values), scale


def minor(frame, prime):
    residues = [[int(v % prime) for v in row] for row in frame.tolist()]
    a = nmod_mat(residues, prime)
    rows, cols = native.pivotcols(a.transpose()), native.pivotcols(a)
    assert len(rows) == len(cols)
    square = nmod_mat([[residues[i][j] for j in cols] for i in rows], prime)
    determinant = int(square.det())
    assert determinant and frame.rank() == len(rows)
    return rows, cols, determinant


def matrix(rows, columns):
    return fmpq_mat([[str(r.get(k, 0)) for k in columns] for r in rows])


def calculate(seed, terms, workers):
    stable = [r for group in native.words_by_shape(10, (2, 2, 2, 2, 2)).values()
              for r in native.singlets(group)]
    assert len(stable) == 370 and terminal_dimension(2) == 162
    raw_basis, _ = albert_basis()
    generators = np.array([g // gcd(*(int(v) for v in g.flat)) for g in raw_basis])
    traces, record, frame = integral_rows(stable, "F4", 4, 162, seed, workers,
                                           basis=generators, sample_terms=terms, return_frame=True)
    recipes = singlets()
    graphs = sorted({g for _, _, r in recipes for g in r})
    ids = {g: i for i, g in enumerate(graphs)}
    coeff = fmpz_mat(len(recipes), len(graphs))
    for i, (_, _, eigenvector) in enumerate(recipes):
        for g, n in eigenvector.items():
            coeff[i, ids[g]] = n
    batches = []
    for start in range(0, frame.P, 2):
        raw, scale = contract(frame, graphs, start, min(start + 2, frame.P))
        batches.append(raw)
    native_values = coeff * fmpz_mat(np.concatenate(batches, axis=1).tolist())
    inactive = sum(int(x) for x, _, _ in recipes)
    assert all(r[0] == (i < inactive) for i, r in enumerate(recipes))
    joined = fmpz_mat(traces.tolist() + native_values.tolist())
    sub = fmpz_mat(joined.tolist()[:len(stable) + inactive])
    print("native trace / inactive / calibration ranks", traces.rank(), sub.rank(), joined.rank(), flush=True)
    assert joined.rank() <= 162
    prime = 1009
    rows, cols, det = minor(sub, prime)
    all_rows, all_cols, all_det = minor(joined, prime)
    basis, index = quartic_pivots()
    five, scores = obstruction_basis(basis, index)
    projected = [full_reduction(reduce_native(r), basis, index, 2) for r in stable]
    keys = sorted(set().union(*(r.keys() for r in projected + five)))
    reductions = matrix(projected + [{} for _ in range(inactive)], keys)
    sampled = fmpq_mat([[sub[i, j] for j in cols] for i in range(sub.nrows())])
    square = fmpq_mat([[sub[i, j] for j in cols] for i in rows])
    selected = fmpq_mat([[reductions[i, j] for j in range(len(keys))] for i in rows])
    residual = reductions - sampled * square.solve(selected)
    relations, rank = residual.rref()
    relation_rows = relations.tolist()[:rank]
    tests = matrix(five, keys)
    intersection = fmpq_mat(relation_rows + (-tests).tolist()).transpose()
    # Compute the rational kernel via an integer numerator and its column basis.
    integer, _ = intersection.numer_denom()
    kernel, n = integer.nullspace()
    coordinates = [[int(kernel[rank+i, j]) for i in range(5)] for j in range(n)]
    coordinates = [[x // gcd(*z) for x in z] for z in coordinates]
    values = [sum(x*y for x, y in zip(z, scores)) for z in coordinates]
    return dict(seed=seed, terms_per_letter=terms, exterior_frame=record,
                defining_singlets=len(stable), formal_Jordan_recipes=len(recipes),
                inactive_Jordan_recipes=inactive, evaluated_graphs=len(graphs),
                albert_inverse_gram_denominator=scale, native_character=162,
                trace_rank=traces.rank(), trace_and_inactive_rank=sub.rank(), calibration_rank=joined.rank(),
                integer_evaluation_injective_on_all_singlets=(joined.rank() == 162),
                test_prime=prime, inactive_minor=dict(rows=rows, columns=cols, determinant=det),
                calibration_minor=dict(rows=all_rows, columns=all_cols, determinant=all_det),
                projected_relation_rank=rank, intersection_dimension=n,
                five_source_coordinates=coordinates, eta_Q1_pairings=list(map(str, values)),
                independent_recipes=[dict(inactive=x, recipe=r) for x, r, _ in recipes])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=736422)
    parser.add_argument("--sample-terms", type=int, default=6)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1 or not 1 <= args.sample_terms <= 52:
        parser.error("enable assertions and choose a positive worker and sample count")
    data = calculate(args.seed, args.sample_terms, args.workers)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print("reduced kernel / intersection / corrected scores", data["projected_relation_rank"],
          data["intersection_dimension"], data["eta_Q1_pairings"])
