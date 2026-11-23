#!/usr/bin/env python3
"""Exact finite seed and AA-transport checks for first Sp relation stability.

The accompanying proof supplies the all-rank argument. This program
certifies the finite seed and checks the literal raising-map transport.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import time

HERE=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,default=HERE.parents[1])
    ap.add_argument('--ranks',type=int,nargs='+',default=[3,4,5,6,7]);args=ap.parse_args()
    base=args.base.resolve()
    spec=importlib.util.spec_from_file_location('baseline_replay',base/'replay.py')
    baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)
    baseline.environment();sys.path.insert(0,str(base/'symplectic/sp4'))
    import paired_relation_basis as p
    from flint import fmpz_mat,nmod_mat
    scalar=(0,0,(0,));pair=(scalar,scalar)
    start=time.monotonic();records=[];previous=None;seed=None
    for rank in args.ranks:
        basis=sorted({pairs for content in p.m.contents((3,3,2*rank-1,0,0),2*rank+2)
                      for pairs in p.paired_content(tuple(sorted(content)))})
        images=[p.raising_pairs(pairs) for pairs in basis]
        target=sorted({v for image in images for v in image});ids={v:i for i,v in enumerate(target)}
        matrix=fmpz_mat(len(target),len(basis))
        for j,row in enumerate(images):
            for word,c in row.items():matrix[ids[word],j]=c
        kernel,dimension=matrix.nullspace()
        transport=None;commutes=None
        if previous is not None and rank==previous[0]+1:
            old_basis,old_images=previous[1:]
            transported=[tuple(sorted(paths+(pair,))) for paths in old_basis]
            transport=sorted(transported)==basis
            lookup=dict(zip(basis,images))
            commutes=all(lookup[new]=={tuple(sorted(paths+(pair,))):c for paths,c in image.items()}
                         for new,image in zip(transported,old_images))
            assert commutes
            if rank>=5:assert transport
        if rank==4:
            reduced,matrix_rank=nmod_mat(matrix.tolist(),1000003).rref()
            pivots=[];j=0
            for i in range(matrix_rank):
                while not reduced[i,j]:j+=1
                pivots.append(j);j+=1
            assert matrix_rank==len(target)==251 and len(basis)==318 and dimension==67
            minor=fmpz_mat([[matrix[i,j] for j in pivots] for i in range(len(target))])
            determinant=int(minor.det());assert determinant
            seed=dict(weight_space_dimension=318,raised_weight_space_dimension=251,
                      raising_rank=251,highest_weight_dimension=67,pivot_columns=pivots,
                      exact_integer_minor_determinant=str(determinant))
        record=dict(rank=rank,pair_weight_space=len(basis),raised_weight_space=len(target),
                    highest_weight_dimension=dimension,AA_divisibility_failures=sum(pair not in x for x in basis),
                    previous_rank_AA_transport_bijective=transport,raising_commutes_with_transport=commutes)
        records.append(record);previous=(rank,basis,images)
    result=dict(passed=True,records=records,rank_four_seed=seed,
                all_rank_proof='See README.md: every nonzero pair monomial at rank>=5 contains AA.',
                scope='first relation degree only; no all-rank Q0 cohomology claim',
                elapsed_seconds=time.monotonic()-start)
    (HERE/'STABILITY_CERTIFICATE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
