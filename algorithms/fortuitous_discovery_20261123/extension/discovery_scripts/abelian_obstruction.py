"""Use Q0 and exact supercommutative trace restriction as upper certificates.

The restriction sorts letters with their original Grassmann signs and
keeps the original trace generators. Odd-length symplectic Cartan traces
vanish. There is no field rescaling and no numerical Cartan sampling.
"""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import argparse
import hashlib
import json
import pickle
import sys
import time

HERE=Path(__file__).resolve().parent
JETS=HERE.parent/'fortuity_discovery_20261123/jets'
sys.path.insert(0,str(JETS))
import pfaffian_relations as pf
from paired_relation_basis import canonical_items
from flint import nmod_mat, fmpq_mat, fmpz_mat
from fractions import Fraction
import math
from fortuity.full.letters import Q_letter, letters_up_to


@lru_cache(None)
def word_restriction(word):
    sign=1
    traces=[]
    for trace in word:
        if len(trace)%2:
            return 0,()
        phase,ordered=canonical_items(trace,pf.is_ferm)
        if not phase:
            return 0,()
        sign*=phase
        traces.append(ordered)
    phase,ordered=canonical_items(traces,lambda tr:sum(map(pf.is_ferm,tr))%2)
    return sign*phase,ordered


def restrict(row):
    result=Counter()
    for word,c in row.items():
        sign,ordered=word_restriction(word)
        if sign:result[ordered]+=sign*c
    return {w:c for w,c in result.items() if c}


def rank_certificate(rows,prime):
    echelon=pf.Echelon(prime);selected=[]
    for i,row in enumerate(rows):
        if echelon.add(row):selected.append(i)
    inverse={j:w for w,j in echelon.words.items()}
    columns=[inverse[j] for j in sorted(echelon.rows)]
    minor=nmod_mat([[rows[i].get(w,0)%prime for w in columns] for i in selected],prime)
    determinant=int(minor.det())
    assert determinant
    return dict(rank=len(selected),rows=selected,columns=columns,determinant_mod_prime=determinant)


def rational_rank_upper(rows,certificate):
    selected=[rows[i] for i in certificate['rows']]
    columns=certificate['columns']
    minor=fmpq_mat([[row.get(w,0) for row in selected] for w in columns])
    rhs=fmpq_mat([[row.get(w,0) for row in rows] for w in columns])
    coefficients=minor.solve(rhs)
    exact_coordinates=[]
    for j,row in enumerate(rows):
        values=[Fraction(str(coefficients[i,j])) for i in range(len(selected))]
        denominator=math.lcm(*(c.denominator for c in values))
        combination=Counter()
        for c,source in zip(values,selected):
            if c:
                multiplier=int(c*denominator)
                for w,n in source.items():combination[w]+=multiplier*n
        assert {w:c for w,c in combination.items() if c}=={w:c*denominator for w,c in row.items()}
        exact_coordinates.append([str(c) for c in values])
    return exact_coordinates


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rank',type=int,choices=[4,5],required=True)
    ap.add_argument('--prime',type=int,default=1000003)
    args=ap.parse_args();rank=args.rank;start=time.monotonic()
    alphabet=letters_up_to(3,3)
    for letter in alphabet:
        normal=Counter()
        for c,left,right in Q_letter(letter):
            sign,pair=canonical_items((left,right),pf.is_ferm)
            if sign:normal[pair]+=c*sign
        assert not any(normal.values()),('Q letter under supercommutation',letter)
    source_path=HERE/f'complete_target_sp{rank}.pkl'
    rows,_=pickle.loads(source_path.read_bytes());assert len(rows)==183
    qrows=[];restrictions=[];combined=[]
    for i,row in enumerate(rows):
        path=HERE/f'outgoing_sp{rank}/row_{i}.pkl'
        digest,image=pickle.loads(path.read_bytes())
        assert digest==hashlib.sha256(pickle.dumps(row,protocol=4)).hexdigest()
        assert not restrict(image),('restriction of Q image',i)
        restriction=restrict(row)
        qrows.append(image);restrictions.append(restriction)
        combined.append({**{('Q0',w):c for w,c in image.items()},
                         **{('Cartan',w):c for w,c in restriction.items()}})
    qcert=rank_certificate(qrows,args.prime)
    ccert=rank_certificate(restrictions,args.prime)
    both=rank_certificate(combined,args.prime)
    integer_minor=fmpz_mat([[combined[i].get(w,0) for w in both['columns']] for i in both['rows']])
    integer_determinant=integer_minor.det()
    assert int(integer_determinant)%args.prime==both['determinant_mod_prime']
    qcoordinates=rational_rank_upper(qrows,qcert)
    ccoordinates=rational_rank_upper(restrictions,ccert)
    result=dict(group=f'Sp({rank})',native_weight=[3,3,2*rank-1,0,0],target_relation_dimension=183,
                prime=args.prime,Q0_rank_lower_bound=qcert['rank'],
                Cartan_rank_lower_bound=ccert['rank'],combined_rank_lower_bound=both['rank'],
                combined_minor_rows=both['rows'],
                combined_minor_columns=[{'map':kind,'word':pf.m.encode(word)} for kind,word in both['columns']],
                combined_minor_determinant_mod_prime=both['determinant_mod_prime'],
                combined_minor_determinant_integer=str(integer_determinant),
                exact_restriction_of_all_target_Q0_images_zero=True,
                exact_letter_restriction_checks=len(alphabet),
                Q0_rational_rank=qcert['rank'],Cartan_rational_rank=ccert['rank'],
                input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                elapsed_seconds=time.monotonic()-start)
    if both['rank']==61:
        result.update(exact_combined_rank=61,common_kernel_dimension=183-61,
                      exact_fortuitous_quotient_dimension=(183-61)-121,
                      target_relation_cohomology_dimension=183-qcert['rank']-121,
                      Cartan_rank_on_Q0_closed_relations=61-qcert['rank'],
                      conditional_on_previously_verified_rational_witness_and_source_rank=True)
    (HERE/f'abelian_obstruction_sp{rank}.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/f'abelian_restrictions_sp{rank}.pkl').write_bytes(pickle.dumps(restrictions,protocol=4))
    (HERE/f'rational_relation_map_coordinates_sp{rank}.json').write_text(json.dumps(dict(
        Q0_basis_rows=qcert['rows'],Q0_coordinates=qcoordinates,
        Cartan_basis_rows=ccert['rows'],Cartan_coordinates=ccoordinates),indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('combined_minor_rows','combined_minor_columns')}),flush=True)


if __name__=='__main__':main()
