"""Extend the two derivative-free F4 classical images using both trace representations.

Each row is the classical variation of a flavour-singlet BMN contraction. Both
representations use identical Lie-algebra coordinates on every evaluation point.
Only independent rows are retained between batches.
"""
import argparse
import hashlib
import random
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
import scalar_neighbors as scalar
from f4_adjoint_traces import adjoint_frame
from f4_bmn_primary import WEIGHT, bmn_convention, translated
from f4_colour import frame
from f4_dressed_jordan import source, variation
from f4_formal_cycles import probe
from f4_joined_boundaries import native
from f4_scalar_correction import evaluate_rows
from f4_singlet_characters import singlets
from mixed33_exact import read


def draw(length, seed):
    rng = random.Random(seed)
    kind = "p" if length == 11 else rng.choice(("f", "pp", "pp"))
    slots = list(range(length))
    rng.shuffle(slots)
    sizes = rng.choice(tuple(scalar.shapes(length)))
    layout = []
    for size in sizes:
        layout.append(tuple(slots[:size]))
        slots = slots[size:]
    return kind, tuple(layout)


def matrices(recipes, dimensions, prime, libs, workers, salt, colours):
    expressions = [translated(scalar.expression(recipe)) for recipe in recipes]
    words, rows = m.flatten(expressions)
    if not words:
        return {}, set(), []
    letters = {letter for w in words for trace in w for letter in trace}
    assert len(letters) == (3 if 0 in dimensions else 6)
    values = {}
    for colour in colours:
        make_frame = frame if colour == 27 else adjoint_frame
        values[colour], frames = probe(rows, words, dimensions, prime, libs[colour],
                                       workers, salt, make_frame)
    return values, letters, frames


def baseline(rec, dimensions, prime, libs, jordan, work, workers, salt):
    n = len(rec["ordinary_recipes"])
    pivots = rec["minor"]["rows"]
    length = rec["length"]
    recipes = [scalar.decode(rec["ordinary_recipes"][i]) for i in pivots if i < n]
    ordinary, alphabet, frames = matrices(recipes, dimensions, prime, libs, workers, salt, [27])
    assert [{k:v for k,v in f.items() if k != "words"} for f in frames] == [
        {k:v for k,v in f.items() if k != "words"} for f in rec["frames"]]
    seeds = [rec["jordan_seeds"][i-n] for i in pivots if i >= n]
    exceptional = [variation(source(0, 4, 12-length, seed), 0) for seed in seeds]
    joined, _ = evaluate_rows(exceptional, alphabet, frames, prime, work, workers, jordan, False)
    result = np.vstack((ordinary[27], joined))
    determinant = int(nmod_mat(result[:,rec["minor"]["columns"]].tolist(), prime).det())
    assert determinant == rec["minor"]["determinant"]
    return result, [("previous", i) for i in pivots]


def run(previous, work, workers, draws):
    work.mkdir(exist_ok=True, parents=True)
    old = read(previous)
    p = old["prime"]
    bmn_convention()
    libs = {n:m.library(work, n, p) for n in (27, 52)}
    jordan = native(work, p, "jordan")
    records = []
    incoming_predecessor = old["records"][1]["minor"]["rank"]
    for rec in old["records"]:
        length = rec["length"]
        salt = 715180 + length*73
        dimensions = singlets("F4", WEIGHT, length+1)
        mat, labels = baseline(rec, dimensions, p, libs, jordan, work, workers, salt)
        limit = min(sum(dimensions.values())-1, 141-incoming_predecessor) if length == 11 else 141-records[0]["minor"]["rank"]
        history = []
        first = 4195819 + length*9301
        for start in range(first, first+draws, 96):
            seeds = list(range(start, min(start+96, first+draws)))
            recipes = [draw(length, seed) for seed in seeds]
            blocks, _, _ = matrices(recipes, dimensions, p, libs, workers, salt, (27,52))
            if not blocks:
                continue
            mat = np.vstack((mat, blocks[27], blocks[52]))
            labels += [(n, seed) for n in (27,52) for seed in seeds]
            minor = m.minor(mat,p)
            chosen = minor["rows"]
            mat = mat[chosen]
            labels = [labels[i] for i in chosen]
            history.append([seeds[-1], len(chosen)])
            assert len(chosen) <= limit
            print("length",length,"rank",len(chosen),"ceiling",limit,flush=True)
            if len(chosen) == limit:
                break
        minor = m.minor(mat,p)
        records.append(dict(length=length,rows=labels,minor=minor,history=history,
                            previous_rank=rec["minor"]["rank"]))
        m.write(work/"trace_images.json",dict(prime=p,input_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
                weight=WEIGHT,dimensions=old["dimensions"],records=records,
                classical_cohomology_upper=141-sum(v["minor"]["rank"] for v in records) if len(records)==2 else None))
    return [r["minor"]["rank"] for r in records]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, default=Path(__file__).with_name("f4")/"primary.json")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--draws", type=int, default=1152)
    args = parser.parse_args()
    assert __debug__ and min(args.workers,args.draws)>0
    print(run(**vars(args)), flush=True)
