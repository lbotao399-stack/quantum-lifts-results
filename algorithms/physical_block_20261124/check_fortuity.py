#!/usr/bin/env python3
"""Exact full stable-to-finite fortuity quotient at the Sp(2) threshold.

The relation space is the entire nullspace of the positive finite Wick form.
Its image is differentiated in the original universal trace algebra.
"""
from compute_metric import *

def exact_rank_rows(rows):
    words=sorted({w for r in rows for w in r});ids={w:i for i,w in enumerate(words)}
    matrix=fmpz_mat(len(rows),len(words))
    for i,row in enumerate(rows):
        denominator=math.lcm(*(Fraction(c).denominator for c in row.values())) if row else 1
        for w,n in row.items():matrix[i,ids[w]]=int(Fraction(n)*denominator)
    print('exact rank input',matrix.nrows(),matrix.ncols(),flush=True)
    rank=matrix.rank()
    print('exact rank result',rank,flush=True)
    return rank

def check_length(length):
    start=time.monotonic()
    centre=read(HERE/f'frames_Y{length}_offset0.json.gz')
    basis=[r for d,rs in sorted(centre['basis'].items(),key=lambda x:int(x[0])) for r in import_rows(rs)]
    derivative=[q(r) for r in basis]
    record=read(HERE/f'finite_metric_Y{length}.json.gz')
    relations=import_matrix(record['null_relations'])
    relrows=[]
    for j in range(relations.ncols()):
        row=Counter()
        for i in range(relations.nrows()):
            c=rational(relations[i,j])
            if c:
                for w,n in derivative[i].items():row[w]+=c*n
        relrows.append({w:n for w,n in row.items() if n})
    stable_rank=exact_rank_rows(derivative)
    relation_rank=exact_rank_rows(relrows)
    harmonic=read(HERE/f'harmonic_Y{length}.json.gz')
    outgoing=import_matrix(harmonic['outgoing_gram'])
    _,finite_rank=outgoing.rref()
    dimension=stable_rank-relation_rank-finite_rank
    incoming=import_matrix(harmonic['incoming_pairing'])
    _,incoming_rank=incoming.rref()
    assert finite_rank+incoming_rank+harmonic['harmonic_dimension']==record['finite_dimension']
    assert 0<=dimension<=harmonic['harmonic_dimension']
    if length==7:assert dimension==1 and harmonic['harmonic_dimension']==1
    result={'passed':True,'length':length,'full_stable_dimension':len(basis),
            'full_finite_dimension':record['finite_dimension'],'all_relation_count':relations.ncols(),
            'stable_outgoing_rank':stable_rank,'relation_image_rank':relation_rank,
            'finite_outgoing_rank':finite_rank,'finite_incoming_rank':incoming_rank,
            'fortuitous_quotient_dimension':dimension,
            'stable_image_in_finite_cohomology_dimension':harmonic['harmonic_dimension']-dimension,
            'full_finite_classical_cohomology_dimension':harmonic['harmonic_dimension'],
            'seconds':time.monotonic()-start,
            'scope':'Complete universal symplectic trace cover, with every finite relation and incoming boundary. No Cartan test is substituted.'}
    print(result,flush=True)
    return result

def main():
    results=[check_length(length) for length in (5,6,7,8)]
    result=dict(results[2]);result['all_native_length_blocks']=results
    result['total_fortuitous_multiplicity']=sum(r['fortuitous_quotient_dimension'] for r in results)
    result['total_stable_image_multiplicity']=sum(r['stable_image_in_finite_cohomology_dimension'] for r in results)
    write(HERE/'fortuity_receipt.json',result)

if __name__=='__main__':main()
