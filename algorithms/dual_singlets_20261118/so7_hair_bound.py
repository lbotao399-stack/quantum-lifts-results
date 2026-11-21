"""Bound the rational boundary rank for the SO(7) level-26 first differential."""
import argparse
from pathlib import Path
import sys

import numpy as np
import mixed_incoming as m
from dual_charge_bands import coordinates
from trace_descent import dimensions, highest_weights

sys.path.append(str(Path(__file__).resolve().parents[1]/"so7_cartan_target_20261118"))
from reduce_incoming import Echelon


def projected_image(weight,length,colours,widths,required,seed,out,workers,output_name,
                    checkpoint=False):
    out.mkdir(exist_ok=True, parents=True)
    group = "Sp" if colours == 4 else "SO"
    previous = dimensions(weight, length+1, colours, group)
    assert set(widths) == set(previous)
    alphabet = sorted({l for c in m.contents(weight,length+1) for l in c})
    formal = highest_weights(weight,length)
    rows = [v for degree in sorted(formal) for v in formal[degree]]
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(rows)).tolist()
    # A minor on any subset of evaluations bounds the rational image rank.
    samples = {}
    for degree,n in widths.items():
        points = (n+degree)//(degree+1)
        samples[degree]=[(min(6,points-j),785941+197*degree+6151*j,
                          list(range(min(6,points-j)*(degree+1)))) for j in range(0,points,6)]
    width = sum(len(cols) for by_degree in samples.values() for _,_,cols in by_degree)
    space = Echelon(width,minors=True)
    lib = m.library(out,colours,m.P)
    seen = 0
    def save_minor(path):
        panels = [dict(rows=p["ids"], columns=list(map(int,p["pivots"])),det=p["determinant"])
                  for p in space.panels]
        assert all(p["det"]%m.P for p in panels)
        assert space.rank == sum(len(p["rows"]) for p in panels)
        result = dict(weight=weight,group=group,colours=colours,prime=m.P,source_length=length,
                      target_dimension=sum(previous.values()),formal_inputs=len(rows),
                      sampled_inputs=seen,permutation_seed=seed,projected_coordinates=width,
                      target_evaluations=samples,closed_rank_lower_bound=space.rank,
                      next_boundary_rank_upper_bound=sum(previous.values())-space.rank,
                      panels=panels)
        m.write(path,result)
        return result
    for i in range(0,len(rows),192):
        indices = order[i:i+192]
        images = [m.q(rows[j]) for j in indices]
        values = np.column_stack([coordinates(images,d,alphabet,s,colours,group,lib,workers)
                                  for d,s in samples.items()])
        space.add(values,indices)
        seen += len(indices)
        if seen%1536 == 0:
            print("closed image rows",seen,space.rank,flush=True)
            if checkpoint:
                save_minor(out/"partial.json.gz")
        if space.rank >= required:
            break
    return save_minor(out/output_name)


def run(out, workers):
    return projected_image((1,1,4,3,3),8,7,{1:4500,3:9000,5:3200},12400,108723,
                           out,workers,"previous.json.gz")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=8)
    args=parser.parse_args()
    run(args.out,args.workers)
