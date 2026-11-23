"""Compare the nine-input F4 inactive primitives with their independent characters."""
import argparse
from fractions import Fraction
import json
from math import gcd
from pathlib import Path

import numpy as np
from flint import fmpz_mat

from jordan_bracket_test import albert_basis, contract, integral_rows, minor, native, singlets
from source_symmetry import centralizer, character, partitions, singlet_trace


def source_dimension(kind):
    ghost, scalar = ((1, 1), (2, 2, 2)) if kind == 2 else ((2, 1), (2, 2, 1))
    dimension = sum(Fraction(character(ghost, x)*character(scalar, y) *
                             singlet_trace(tuple(sorted((*x, *y, 1)))),
                             centralizer(x)*centralizer(y))
                    for x in partitions(sum(ghost)) for y in partitions(sum(scalar)))
    assert dimension.denominator == 1
    return int(dimension)


def sector_of(row):
    word = next(iter(row))
    (second,) = [letter for trace in word for letter in trace
                 if letter[0]+letter[1]+len(letter[2])==2]
    return sum(second[:2])


def calculate(kind, seed, workers):
    upper = source_dimension(kind)
    shapes = native.words_by_shape(9, (2, 2, 2, 2, 2))
    stable = [r for group in shapes.values() if sector_of({group[0]:1}) == kind
              for r in native.singlets(group)]
    raw_basis, _ = albert_basis()
    generators = np.array([g // gcd(*(int(v) for v in g.flat)) for g in raw_basis])
    traces, record, frame = integral_rows(stable, "F4", 3, upper, seed, workers,
                                           basis=generators, sample_terms=6, return_frame=True)
    rows = [(data, row) for _, data, row in singlets(sources=True)
            if sum(i in data[1] for i in data[-1]) == kind]
    graphs = sorted({g for _, row in rows for g in row})
    ids = {g: i for i, g in enumerate(graphs)}
    coeff = fmpz_mat(len(rows), len(graphs))
    for i, (_, row) in enumerate(rows):
        for g, n in row.items():
            coeff[i, ids[g]] = n
    joined = traces
    if rows:
        batches = [contract(frame, graphs, i, min(i+2, frame.P))[0] for i in range(0, frame.P, 2)]
        extra = coeff * fmpz_mat(np.concatenate(batches, axis=1).tolist())
        joined = fmpz_mat(traces.tolist() + extra.tolist())
    selected, cols, determinant = minor(joined, 1009)
    assert len(selected) <= upper
    return dict(dotted_degree_on_second_jet=kind,native_character=upper,
                formal_trace_singlets=len(stable), evaluated_trace_rank=traces.rank(),
                formal_Jordan_singlets=len(rows),distinct_graphs=len(graphs),
                combined_evaluation_rank=len(selected),external_frame=record,
                minor_prime=1009,minor_determinant=determinant,
                minor_rows=selected,minor_columns=cols,Jordan_recipes=[r for r,_ in rows])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind",type=int,choices=(1,2),required=True)
    parser.add_argument("--seed",type=int,default=921554)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("enable assertions and use a positive worker count")
    data=calculate(args.kind,args.seed,args.workers)
    args.out.write_text(json.dumps(data,indent=2)+"\n")
    print(args.kind, data["evaluated_trace_rank"],data["combined_evaluation_rank"],data["native_character"])
