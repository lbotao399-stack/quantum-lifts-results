"""Rank a BMN classical differential by evaluating its image without framing the codomain."""
import argparse
from pathlib import Path

import numpy as np

import mixed_incoming as m
from dual_charge_bands import character_frame, coordinates
from orthosymplectic_classical_page import minor
from trace_descent import dimensions, highest_weights


def image_samples(columns, degree, seed):
    points = (columns + degree) // (degree + 1)
    samples = []
    for start in range(0, points, 6):
        count = min(6, points - start)
        samples.append((count, seed + 6151 * start, list(range(count * (degree+1)))))
    return samples


def run(weight, length, group, workers, seed, out):
    out.mkdir(parents=True, exist_ok=True)
    colours = {"Sp": 4, "SO": 5}[group]
    w = (0, 0, *weight)
    source_dims = dimensions(w, length, colours, group)
    target_dims = dimensions(w, length+1, colours, group)
    source_alphabet, target_alphabet = [
        sorted({l for content in m.contents(w,y) for l in content})
        for y in (length, length+1)]
    lib = m.library(out, colours, m.P)
    formal = highest_weights(w,length)
    record, variations = {}, []
    for degree, n in source_dims.items():
        _, data = character_frame(formal[degree], n, degree, source_alphabet,
                                   colours, group, lib, workers, seed+317*degree)
        basis = [formal[degree][j] for j in data["minor_rows"]]
        variations.extend(m.q(row) for row in basis)
        record[degree] = data
        print("source", degree, len(basis), len(formal[degree]), flush=True)
    # Each variation of a BMN trace changes the fermion count by one.
    # Only a lower bound is inferred from sampling any codomain coordinates.
    target_records, matrices = {}, []
    for degree,n in target_dims.items():
        supported = sum(any(m.nferm(word) == degree for word in row) for row in variations)
        if not supported:
            continue
        samples = image_samples(min(n,supported)+12, degree, seed+80113+593*degree)
        target_records[degree] = {"character":n, "supported_rows":supported, "samples":samples}
        parts = [coordinates(variations[i:i+384], degree, target_alphabet, samples,
                             colours, group, lib, workers) for i in range(0,len(variations),384)]
        matrices.append(np.vstack(parts))
        print("image", degree, matrices[-1].shape, flush=True)
    image = np.column_stack(matrices)
    result = dict(group=group, colours=colours, weight=weight, length=length,
                  prime=m.P, seed=seed, source_dimension=sum(source_dims.values()),
                  source_frames=record, target_evaluations=target_records,
                  rank_minor=minor(image,m.P))
    m.write(out/"image.json.gz", result)
    print("outgoing rank at least",result["rank_minor"]["rank"],flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", nargs=3, type=int, default=(10,4,4))
    parser.add_argument("--length", type=int, default=10)
    parser.add_argument("--group", choices=("Sp","SO"), default="Sp")
    parser.add_argument("--workers",type=int, default=4)
    parser.add_argument("--seed",type=int, default=420569)
    parser.add_argument("--out",type=Path, required=True)
    args=parser.parse_args()
    assert __debug__
    run(tuple(args.weight), args.length, args.group, args.workers, args.seed, args.out)
