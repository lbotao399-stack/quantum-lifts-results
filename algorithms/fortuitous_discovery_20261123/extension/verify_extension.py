#!/usr/bin/env python3
"""Exact upper certificates for the two previously certified Sp classes.

Install this directory as `extension/` inside the four-certificate package,
or pass --base explicitly. No original discovery workspace is required.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import pickle
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent


def tuples(value):
    return tuple(tuples(v) for v in value) if isinstance(value,list) else value


def normal_order(items,odd):
    fermions=[v for v in items if odd(v)]
    if len(set(fermions))!=len(fermions):return 0,()
    phase=(-1)**sum(a>b for i,a in enumerate(fermions) for b in fermions[i+1:])
    return phase,tuple(sorted(items))


def semantic_digest(row,encode):
    data=[[str(row[word]),encode(word)] for word in sorted(row)]
    return hashlib.sha256(json.dumps(data,separators=(',',':')).encode()).hexdigest()


def verify_coordinates(rows,basis_indices,coordinates,fmpq_mat):
    basis=[rows[i] for i in basis_indices]
    assert len(coordinates)==len(rows)
    for row,entries in zip(rows,coordinates):
        coefficients=list(map(Fraction,entries))
        assert len(coefficients)==len(basis)
        denominator=math.lcm(*(c.denominator for c in coefficients))
        combined=Counter()
        for coefficient,source in zip(coefficients,basis):
            if coefficient:
                multiplier=int(coefficient*denominator)
                for word,c in source.items():combined[word]+=multiplier*c
        assert {w:c for w,c in combined.items() if c}=={w:c*denominator for w,c in row.items()}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',type=Path,default=HERE.parent)
    parser.add_argument('--rank',type=int,choices=[4,5],required=True)
    parser.add_argument('--upper-only',action='store_true',help='Reuse the separately verified baseline witness.')
    args=parser.parse_args();base=args.base.resolve();rank=args.rank
    start=time.monotonic()
    spec=importlib.util.spec_from_file_location('baseline_replay',base/'replay.py')
    baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)
    baseline.environment()
    folder=baseline.hydrate(rank);sys.path.insert(0,str(folder))
    import pfaffian_relations as pf
    import path_relation_completion as paths
    from fortuity.full.letters import Q_letter,letters_up_to
    from flint import fmpq_mat,fmpz_mat
    if not args.upper_only:
        baseline.symplectic(rank)

    certificate=json.loads((HERE/f'certificate_sp{rank}.json').read_text())
    coordinates=json.loads((HERE/f'rational_map_coordinates_sp{rank}.json').read_text())
    additions=json.loads((HERE/f'additional_identities_sp{rank}.json').read_text())
    old_certificate=json.loads((folder/f'sp{rank}_discovery_certificate.json').read_text())
    assert old_certificate['exact_rational_source_relation_image_rank']==121
    assert old_certificate['finite_rational_cohomology_lower_bound']==1
    name='target_snapshot_181.pkl' if rank==4 else 'target_sp5_snapshot_181.pkl'
    rows,_=pickle.loads((folder/name).read_bytes());assert len(rows)==181
    expand=lru_cache(maxsize=16)(paths.expand_paths)
    for record in additions['identities']:
        polynomial=Counter()
        for coefficient,graph in record['terms']:
            graph=tuples(graph)
            if record['kind']=='paths':
                image=expand(graph,rank)
            else:
                assert record['kind']=='paths_cycle'
                graph_paths,cycle=graph
                image=Counter()
                for word,c in expand(graph_paths,rank).items():
                    sign,new_word=pf.m.so_canon(word+(cycle,))
                    if sign:image[new_word]+=sign*c
            for word,c in image.items():polynomial[word]+=coefficient*c
        row=pf.primitive({w:c for w,c in polynomial.items() if c})
        assert semantic_digest(row,pf.m.encode)==record['polynomial_sha256']
        assert not pf.spin_action(row,True)
        rows.append(row)
    assert len(rows)==183
    echelon=pf.Echelon()
    for row in rows:assert echelon.add(row)
    dimensions,_=pf.profile_inventory(rank,(3,3,2*rank-1,0,0),2*rank+4)
    assert sum(dimensions.values())==183

    alphabet=letters_up_to(3,3)
    for letter in alphabet:
        normal=Counter()
        for c,left,right in Q_letter(letter):
            sign,pair=normal_order((left,right),pf.is_ferm)
            if sign:normal[pair]+=sign*c
        assert not any(normal.values())

    @lru_cache(None)
    def abelian_word(word):
        sign=1;traces=[]
        for trace in word:
            if len(trace)%2:return 0,()
            phase,ordered=normal_order(trace,pf.is_ferm)
            if not phase:return 0,()
            sign*=phase;traces.append(ordered)
        phase,ordered=normal_order(traces,lambda tr:sum(map(pf.is_ferm,tr))%2)
        return sign*phase,ordered

    def restriction(row):
        result=Counter()
        for word,c in row.items():
            sign,ordered=abelian_word(word)
            if sign:result[ordered]+=sign*c
        return {w:c for w,c in result.items() if c}

    qrows=[pf.m.q(row) for row in rows]
    cartan=[restriction(row) for row in rows]
    assert all(not restriction(row) for row in qrows)
    verify_coordinates(qrows,coordinates['Q0_basis_rows'],coordinates['Q0_coordinates'],fmpq_mat)
    verify_coordinates(cartan,coordinates['Cartan_basis_rows'],coordinates['Cartan_coordinates'],fmpq_mat)
    # The integer minor gives the lower bound; the rational coordinates
    # above give both individual upper bounds. The baseline 122 closed
    # vectors give the combined upper bound.
    selected=certificate['combined_minor_rows']
    columns=[(entry['map'],pf.m.decode(entry['word'])) for entry in certificate['combined_minor_columns']]
    minor=fmpz_mat([[(qrows if kind=='Q0' else cartan)[i].get(word,0)
                     for kind,word in columns] for i in selected])
    determinant=int(minor.det())
    assert minor.nrows()==minor.ncols()==61
    assert determinant==int(certificate['combined_minor_determinant_integer'])!=0
    factorization=certificate['combined_minor_determinant_factorization']
    assert determinant==factorization['sign']*math.prod(p**e for p,e in factorization['factors'])
    # Check the two lower rank witnesses as well as their rational spans.
    for image_rows,basis_indices,expected in [(qrows,coordinates['Q0_basis_rows'],45),
                                             (cartan,coordinates['Cartan_basis_rows'],24)]:
        e=pf.Echelon()
        for i in basis_indices:assert e.add(image_rows[i])
        assert len(e.rows)==expected
    result=dict(passed=True,group=f'Sp({rank})',native_weight=[3,3,2*rank-1,0,0],
                exact_target_identity_dimension=183,Q0_rank=45,stable_abelianization_rank=24,
                combined_rank=61,integer_minor_determinant=str(determinant),
                relation_cohomology_dimension=183-45-121,
                stable_abelianization_rank_on_relation_cohomology=61-45,
                fortuitous_dimension=(183-61)-121,
                full_finite_sector_cohomology_dimension=1,
                full_finite_dimension_uses_previously_certified_stable_centre_acyclicity=True,
                native_letter_Q0_abelianization_checks=len(alphabet),
                baseline_witness_replayed=not args.upper_only,
                arithmetic='Q and exact integer determinants',classical_Q0_only=True,
                elapsed_seconds=time.monotonic()-start)
    output=HERE/'replay_outputs';output.mkdir(exist_ok=True)
    (output/f'sp{rank}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
