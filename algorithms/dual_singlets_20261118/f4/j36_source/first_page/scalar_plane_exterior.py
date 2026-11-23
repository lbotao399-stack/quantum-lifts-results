"""Prolong the one-ghost Ward quotient to two equal dotted ghosts on a scalar plane."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
from itertools import combinations_with_replacement
import json
from pathlib import Path
from random import Random

from flint import nmod_mat
import numpy as np

from scalar_plane_ward import blocks, coordinate_groups, plane, tableaux


def load_record(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw)


def one_kernel(data, prime):
    d, m, seed, power = (data[k] for k in ("scalar_colour_dimension", "grading_modulus",
                                             "seed", "determinant_degree"))
    _, degrees, grades, _ = plane(d, m, seed)
    _, column_groups = coordinate_groups(d, power, m, degrees, grades)
    _, scalars = tableaux(d, power)
    position = {key: i for i, key in enumerate(scalars)}
    nullity = data["columns"] - data["rank_lower_bound"]
    result = np.zeros((len(scalars), 52, nullity), dtype=np.uint16)
    kernel_grades = []
    ranks = []
    for record, (w, height, width, entries) in zip(data["blocks"], blocks(d, power, m, seed), strict=True):
        assert (w, height, width) == (record["weight"], record["rows"], record["columns"])
        rows = record["row_pivots"] if record["row_pivots"] is not None else list(range(height))
        cols = record["column_pivots"] if record["column_pivots"] is not None else list(range(width))
        irow = {v: i for i, v in enumerate(rows)}
        icol = {v: i for i, v in enumerate(cols)}
        free = [j for j in range(width) if j not in icol]
        ifree = {j: i for i, j in enumerate(free)}
        n, f = len(rows), len(free)
        assert n == len(cols)
        a, b = nmod_mat(n, n, prime), nmod_mat(n, f, prime)
        for i, j, x in entries:
            if i not in irow:
                continue
            r = irow[i]
            if j in icol:
                c = icol[j]
                a[r, c] = a[r, c] + x
            else:
                c = ifree[j]
                b[r, c] = b[r, c] - x
        assert int(a.det()) == record["minor_determinants"][data["primes"].index(prime)] != 0
        if not f:
            ranks.append(n)
            continue
        solved = a.solve(b)
        basis = np.zeros((width, f), dtype=np.int64)
        basis[cols] = np.asarray(solved.tolist(), dtype=np.int64)
        basis[free] = np.eye(f, dtype=np.int64)
        # Check every original equation, including those outside the minor.
        residual = np.zeros((height, f), dtype=np.int64)
        for i, j, value in entries:
            residual[i] += value * basis[j]
        assert not np.any(residual % prime)
        start = len(kernel_grades)
        for j, (alpha, u) in enumerate(column_groups[w]):
            result[position[alpha], u, start:start+f] = basis[j]
        kernel_grades.extend([w] * f)
        ranks.append(n)
    assert len(kernel_grades) == nullity
    scalar_grades = [sum(degrees[i] for i in a) % m for a in scalars]
    return result, kernel_grades, grades, scalar_grades, ranks


def symmetric_test(job):
    total, m, prime, frame, kernel_grades, grades, scalar_grades = job
    unknowns = [(a, v) for a, w in enumerate(kernel_grades) for v in range(52)
                if (w + grades[v]) % m == total]
    index = {key: i for i, key in enumerate(unknowns)}
    order = list(range(len(scalar_grades)))
    Random(42836).shuffle(order)
    width = len(unknowns)
    # Store a basis of solutions after each group of scalar coordinates.
    remaining = None
    dimensions = [width]
    stops = []
    for start in range(0, len(order), 4):
        selected = order[start:start+4]
        rows = []
        for s in selected:
            by_colour = [np.flatnonzero(v).tolist() for v in frame[s]]
            for u, v in combinations_with_replacement(range(52), 2):
                if (scalar_grades[s] + grades[u] + grades[v]) % m != total:
                    continue
                row = [(index[a, v], int(frame[s, u, a])) for a in by_colour[u]]
                if u != v:
                    row.extend((index[a, u], int(frame[s, v, a])) for a in by_colour[v])
                if row:
                    rows.append(row)
        if rows:
            matrix = nmod_mat(len(rows), width, prime)
            for i, row in enumerate(rows):
                for j, v in row:
                    matrix[i, j] = v
            reduced = matrix if remaining is None else matrix * remaining
            new_null, dim = reduced.nullspace()
            dimensions.append(dim)
            stops.append(start + len(selected))
            if dim == 0:
                break
            basis = nmod_mat(reduced.ncols(), dim, prime)
            for i in range(reduced.ncols()):
                for j in range(dim):
                    basis[i, j] = new_null[i, j]
            assert reduced * basis == nmod_mat(reduced.nrows(), dim, prime)
            remaining = basis if remaining is None else remaining * basis
    return dict(weight=total, initial=width, tested_scalar_coordinates=stops,
                remaining_dimensions=dimensions, remaining=dimensions[-1])


def calculate(data, prime, workers):
    frame, kernel_grades, grades, scalar_grades, ranks = one_kernel(data, prime)
    m = data["grading_modulus"]
    jobs = [(w, m, prime, frame, kernel_grades, grades, scalar_grades) for w in range(m)]
    if workers == 1:
        results = list(map(symmetric_test, jobs))
    else:
        with ProcessPoolExecutor(workers) as pool:
            results = list(pool.map(symmetric_test, jobs))
    return dict(prime=prime, scalar_coordinates=len(scalar_grades), one_ghost_block_ranks=ranks,
                one_ghost_dual_dimension=len(kernel_grades),
                alternating_two_ghost_dual_dimension=sum(r["remaining"] for r in results),
                blocks=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("assertions and a positive worker count are required")
    result = calculate(load_record(args.record), args.prime, args.workers)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(result["prime"], result["one_ghost_dual_dimension"],
          result["alternating_two_ghost_dual_dimension"])
