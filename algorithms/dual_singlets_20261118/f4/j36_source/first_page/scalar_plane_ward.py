"""Test the scalar-ghost relations on a fixed colour plane in the F4 relative complex."""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import combinations
import json
from pathlib import Path
from random import Random

from flint import fmpz_mat, nmod_mat
import numpy as np

from flag_boundaries import adjoint
from six_scalar import straightening


def plane(size, modulus, seed):
    actions = adjoint()
    weights = np.diagonal(actions[:4], axis1=1, axis2=2).T
    rng = Random(seed)
    if modulus == 1:
        grades = [0] * 52
        degrees = [0] * size
    else:
        grading = [rng.randrange(1, modulus) for _ in range(4)]
        grades = [int(w @ grading) % modulus for w in weights]
        pool = sorted(set(grades))
        degrees = rng.sample(pool, size) if size <= len(pool) else [
            pool[i % len(pool)] for i in range(size)]
    coefficients = [[rng.randint(-4, 4) if g == d else 0 for g in grades] for d in degrees]
    assert fmpz_mat(coefficients).rank() == size
    matrices = np.array(coefficients, dtype=np.int64) @ actions.reshape(52, -1)
    brackets = []
    for degree, a in zip(degrees, matrices.reshape(size, 52, 52)):
        cols = [[] for _ in range(52)]
        for k, u in zip(*np.nonzero(a)):
            assert (grades[u] + degree - grades[k]) % modulus == 0
            cols[int(k)].append((int(u), int(a[k, u])))
        brackets.append(cols)
    return coefficients, degrees, grades, brackets


@cache
def straighten(left, right):
    if not left:
        return ((right, 1),)
    if all(i <= j for i, j in zip(left, right)):
        return ((left + right, 1),)
    if all(j <= i for i, j in zip(left, right)):
        return ((right + left, 1),)
    alphabet = sorted(set(left + right))
    relative = {x: i for i, x in enumerate(alphabet)}
    a, b = (tuple(relative[i] for i in col) for col in (left, right))
    return tuple((tuple(alphabet[i] for i in cols), n) for cols, n in straightening(a, b))


def tableaux(size, power):
    triples = list(combinations(range(size), 3))
    # The pair in a source is completed by the marked row [U,X].
    pairs = list(combinations(range(size), 2))
    spectators = [()] if power == 1 else triples
    source = [(a, b) for a in spectators for b in pairs if all(i <= j for i, j in zip(a, b))]
    target = [a + b for a in spectators for b in triples if all(i <= j for i, j in zip(a, b))]
    return source, target


def coordinate_groups(size, power, modulus, degrees, grades):
    source, target = tableaux(size, power)
    rows, cols = defaultdict(list), defaultdict(list)
    def grade(indices, colour):
        return (sum(degrees[i] for i in indices) + grades[colour]) % modulus
    for a, b in source:
        for k in range(52):
            rows[grade(a + b, k)].append((a, b, k))
    for t in target:
        for u in range(52):
            cols[grade(t, u)].append((t, u))
    assert sum(map(len, rows.values())) == 52 * len(source)
    assert sum(map(len, cols.values())) == 52 * len(target)
    return rows, cols


def blocks(size, power, modulus, seed):
    _, degrees, grades, brackets = plane(size, modulus, seed)
    rows, cols = coordinate_groups(size, power, modulus, degrees, grades)
    for weight in sorted(rows.keys() | cols.keys()):
        indices = {key: i for i, key in enumerate(cols[weight])}
        entries = []
        for row, (left, pair, k) in enumerate(rows[weight]):
            for j in range(size):
                if j in pair or not brackets[j][k]:
                    continue
                ordered = (*pair, j)
                phase = (-1) ** sum(x > j for x in pair)
                for t, value in straighten(left, tuple(sorted(ordered))):
                    for u, v in brackets[j][k]:
                        entries.append((row, indices[t, u], phase * value * v))
        yield weight, len(rows[weight]), len(cols[weight]), entries


def pivot_columns(matrix):
    reduced, rank = matrix.rref()
    cols = []
    j = 0
    for i in range(rank):
        while not reduced[i, j]:
            j += 1
        cols.append(j)
        j += 1
    return cols


def minor(record):
    weight, height, width, entries, primes = record
    first = nmod_mat(height, width, primes[0])
    for i, j, x in entries:
        first[i, j] = first[i, j] + x
    rows = pivot_columns(first.transpose())
    if len(rows) == width:
        cols = list(range(width))
    else:
        selected = {v: i for i, v in enumerate(rows)}
        short = nmod_mat(len(rows), width, primes[0])
        for i, j, x in entries:
            if i in selected:
                short[selected[i], j] = short[selected[i], j] + x
        cols = pivot_columns(short)
    assert len(rows) == len(cols)
    irow = {v: i for i, v in enumerate(rows)}
    icol = {v: i for i, v in enumerate(cols)}
    determinants = []
    for p in primes:
        square = nmod_mat(len(rows), len(cols), p)
        for i, j, x in entries:
            if i in irow and j in icol:
                a, b = irow[i], icol[j]
                square[a, b] = square[a, b] + x
        determinants.append(int(square.det()))
    assert determinants[0]
    # An omitted pivot list means all coordinates in their natural order.
    return dict(weight=weight, rows=height, columns=width, rank_lower_bound=len(rows),
                row_pivots=None if len(rows) == height else rows,
                column_pivots=None if len(cols) == width else cols,
                minor_determinants=determinants)


def calculate(size, power, modulus, seed, workers, primes):
    coefficients, degrees, _, _ = plane(size, modulus, seed)
    jobs = ((*b, primes) for b in blocks(size, power, modulus, seed))
    if workers == 1:
        results = list(map(minor, jobs))
    else:
        with ProcessPoolExecutor(workers) as pool:
            results = list(pool.map(minor, jobs))
    return dict(scalar_colour_dimension=size, determinant_degree=power, seed=seed,
                grading_modulus=modulus, plane_grades=degrees, coefficients=coefficients,
                primes=primes, rows=sum(x["rows"] for x in results),
                columns=sum(x["columns"] for x in results),
                rank_lower_bound=sum(x["rank_lower_bound"] for x in results), blocks=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--power", type=int, choices=(1, 2), default=1)
    parser.add_argument("--grading", type=int, default=1)
    parser.add_argument("--seed", type=int, default=361144)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1 or args.size < 3 or args.grading < 1:
        parser.error("assertions and positive dimensions are required")
    result = calculate(args.size, args.power, args.grading, args.seed, args.workers, [1009, 1013])
    if args.out:
        args.out.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    print(json.dumps({k: result[k] for k in ("scalar_colour_dimension", "determinant_degree",
                                              "rows", "columns", "rank_lower_bound")}))
