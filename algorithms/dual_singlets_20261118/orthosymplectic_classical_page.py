"""Compute adjacent classical minors and, for a one-dimensional quotient, its two witnesses."""
import argparse
import hashlib
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from dual_charge_bands import character_frame, coordinates
from filtered_rank import pivots
from trace_descent import dimensions, highest_weights


def frame(weight, length, group, colours, seed, library, workers):
    sizes = dimensions(weight, length, colours, group)
    rows = highest_weights(weight, length)
    alphabet = sorted({letter for content in m.contents(weight, length) for letter in content})
    basis, samples, records = [], {}, {}
    for degree, size in sizes.items():
        samples[degree], records[degree] = character_frame(
            rows[degree], size, degree, alphabet, colours, group, library, workers,
            seed + 317 * degree)
        basis += [rows[degree][i] for i in records[degree]["minor_rows"]]
    return basis, samples, records, alphabet


def matrix(rows, samples, alphabet, group, colours, library, workers, variation=False):
    width = sum(len(cols) for values in samples.values() for _, _, cols in values)
    result = np.empty((len(rows), width), dtype=np.int64)
    for start in range(0, len(rows), 384):
        batch = rows[start:start+384]
        if variation:
            batch = [m.q(row) for row in batch]
        if width:
            result[start:start+len(batch)] = np.column_stack([
                coordinates(batch, degree, alphabet, probes, colours, group, library, workers)
                for degree, probes in samples.items()])
    return result


def minor(matrix, prime):
    columns = pivots(matrix, prime)
    if not columns:
        return dict(rank=0, rows=[], columns=[], determinant=1)
    rows = pivots(matrix[:, columns].T, prime)
    square = matrix[np.ix_(rows, columns)]
    determinant = int(nmod_mat(square.tolist(), prime).det())
    assert determinant and len(rows) == len(columns)
    return dict(rank=len(rows), rows=rows, columns=columns, determinant=determinant,
                minor_sha256=hashlib.sha256(square.astype("<u4").tobytes()).hexdigest())


def one_class(cochains, incoming, outgoing, before, after, prime):
    # Images of the outgoing independent rows give a complement to the cycles.
    stacked = np.vstack((incoming[before["rows"]], cochains[after["rows"]]))
    reduced, rank = nmod_mat(stacked.tolist(), prime).rref()
    assert rank == len(cochains) - 1
    cols = []
    col = 0
    for i in range(rank):
        while not reduced[i, col]:
            col += 1
        cols.append(col)
        col += 1
    free, = set(range(len(cochains))).difference(cols)
    functional = np.zeros(len(cochains), dtype=np.int64)
    functional[free] = 1
    functional[cols] = [(-int(reduced[i, free])) % prime for i in range(rank)]
    assert len(cochains)*(prime-1)**2 < 2**63
    assert not np.any(incoming @ functional % prime)
    values = cochains @ functional % prime
    i = int(np.flatnonzero(values)[0])
    cycle = np.zeros(len(cochains), dtype=np.int64)
    cycle[i] = 1
    if after["rank"]:
        square = nmod_mat(outgoing[np.ix_(after["rows"], after["columns"])].T.tolist(), prime)
        rhs = nmod_mat([[int(v)] for v in outgoing[i, after["columns"]]], prime)
        cycle[after["rows"]] = [(-int(v)) % prime for v in square.solve(rhs).entries()]
    assert not np.any(cycle @ outgoing % prime)
    value = int(cycle @ values % prime)
    assert value
    functional = functional * pow(value, -1, prime) % prime
    assert int((cycle @ cochains % prime) @ functional % prime) == 1
    return dict(cycle_basis_coefficients=cycle.tolist(), target_functional=functional.tolist(),
                pairing=1, classical_closure_nonzero=0, paired_boundary_rows=len(incoming))


def run(weight, length, group, colours, seed, out, workers):
    out.mkdir(parents=True, exist_ok=True)
    library = m.library(out, colours, m.P)
    centre, probes, data, alphabet = frame(weight, length, group, colours, seed, library, workers)
    _, following, following_data, next_letters = frame(
        weight, length+1, group, colours, seed+19609, library, workers)
    previous = highest_weights(weight, length-1)
    parents = [r for d in sorted(previous) for r in previous[d]]
    incoming = matrix(parents, probes, alphabet, group, colours, library, workers, True)
    outgoing = matrix(centre, following, next_letters, group, colours, library, workers, True)
    before, after = minor(incoming, m.P), minor(outgoing, m.P)
    quotient = len(centre)-before["rank"]-after["rank"]
    assert quotient >= 0
    result = dict(weight=weight, length=length, group=group, colours=colours, prime=m.P,
                  seed=seed, centre_frames=data, following_frames=following_data,
                  formal_parents=len(parents), cochains=len(centre), incoming=before, outgoing=after,
                  modular_cohomology=quotient, rational_upper_bound=quotient)
    if quotient == 1:
        basis_values = matrix(centre, probes, alphabet, group, colours, library, workers)
        result["representative"] = one_class(basis_values, incoming, outgoing, before, after, m.P)
    m.write(out/"classical_page.json.gz", result)
    print(group, weight, length, len(centre), before["rank"], after["rank"], quotient, flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, nargs=5, required=True)
    parser.add_argument("--length", type=int, required=True)
    parser.add_argument("--group", choices=("SO", "Sp"), required=True)
    parser.add_argument("--colours", type=int, required=True)
    parser.add_argument("--seed", type=int, default=812329)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert __debug__ and ((args.group=="SO" and args.colours%2) or
                          (args.group=="Sp" and args.colours%2==0))
    run(tuple(args.weight), args.length, args.group, args.colours, args.seed, args.out, args.workers)
