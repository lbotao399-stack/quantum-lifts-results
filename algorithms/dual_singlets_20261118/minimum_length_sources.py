"""Test classical sources in the shortest possible length at each derivative grade."""
import argparse
import math
from pathlib import Path

import numpy as np
from flint import nmod_mat

import trace_descent as descent
from mixed_incoming import P, nferm, pivotcols, write
from terminal_interval_ranks import sample


def source_bound(flavours, spin, grade, size, group, seed):
    weight = (grade+spin, grade) + tuple(x-grade for x in flavours)
    length = weight[2]
    dimensions = descent.dimensions(weight, length, size, group)
    assert not descent.dimensions(weight, length-1, size, group)
    dimension = sum(dimensions.values())
    parity = (1 if group=="SU+" else -1) if group.startswith("SU") else None
    formal = descent.highest_weights(weight, length, parity)
    rows = [row for degree in sorted(formal) for row in formal[degree]]
    assert len(rows) >= dimension
    images = [descent.action(row,0,size,group) for row in rows]
    by_degree = {}
    for word in sorted({w for row in images for w in row}):
        by_degree.setdefault(nferm(word),[]).append(word)
    for row in rows:
        assert all(0 in letter[2] for word in row for trace in word for letter in trace)
    restricted = {degree:[{w:n for w,n in row.items() if nferm(w)==degree} for row in images]
                  for degree in by_degree}
    frame = np.empty((len(rows),0),dtype=np.int64)
    probes = []
    piv = []
    previous = 0
    rank = 0
    stagnant = 0
    for batch in range(max(24,dimension+24)):
        if rank==dimension or not by_degree:
            break
        # Sample all fermion sectors on each pass, retaining their common
        # source coordinates. Every column is evaluation at an integer point.
        points = min(8,max(2,math.ceil((dimension-rank)/sum(k+1 for k in by_degree))))
        for degree, words in by_degree.items():
            point_seed = seed+6151*batch+191*degree
            values = sample(words,degree,size,points,point_seed,group)
            block = descent.multiply_rows(restricted[degree],words,values)
            frame = np.column_stack((frame,np.asarray(block%P,dtype=np.int64)))
            probes.append([degree,points,point_seed])
        if frame.shape[1]<dimension and batch%5!=4:
            continue
        matrix = nmod_mat(frame.tolist(),P)
        piv = pivotcols(matrix)
        rank = len(piv)
        assert rank<=dimension
        stagnant = stagnant+1 if rank==previous else 0
        previous = rank
        if stagnant >= 2:
            break
    if frame.shape[1] and rank<dimension:
        piv=pivotcols(nmod_mat(frame.tolist(),P))
        rank=len(piv)
    selected_rows = []
    determinant = 1
    if rank:
        selected_rows=pivotcols(nmod_mat(frame[:,piv].T.tolist(),P))
        minor=frame[np.ix_(selected_rows,piv)]
        determinant=int(nmod_mat(minor.tolist(),P).det())
        assert determinant
    return dict(group=group,colours=size,flavours=flavours,spin=spin,grade=grade,
                weight=weight,source_length=length,level=2*sum(flavours)+3*spin,
                source_charges=descent.m.sector(weight,length),
                source_dimension=dimension,preceding_dimension=0,formal_sources=len(rows),
                prime=P,probes=probes,minor_rows=selected_rows,minor_columns=piv,
                determinant_modulo_prime=determinant,classical_rank_lower=rank,
                rational_cohomology_upper=dimension-rank)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flavours",nargs=3,type=int,required=True)
    p.add_argument("--spin",type=int,default=0)
    p.add_argument("--grade",type=int,required=True)
    p.add_argument("--group",choices=["SO","Sp","SU+","SU-"],required=True)
    p.add_argument("--size",type=int,required=True)
    p.add_argument("--seed",type=int,default=81173)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert __debug__ and a.flavours==sorted(a.flavours,reverse=True) and 0<=a.grade<=min(a.flavours)
    record=source_bound(a.flavours,a.spin,a.grade,a.size,a.group,a.seed)
    write(a.out,record)
    print(record["source_dimension"],record["classical_rank_lower"],record["rational_cohomology_upper"])
