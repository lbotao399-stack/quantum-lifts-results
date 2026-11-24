#!/usr/bin/env python3
"""Compute complete radial forms, finite null relations, and harmonic states."""
from wick import *
from flint import nmod_mat, fmpq, fmpq_mat
import argparse

def rational(value):return Fraction(str(value))
def zero(matrix):return matrix==fmpq_mat(matrix.nrows(),matrix.ncols())

def rational_nullspace(matrix):
    rr,rank=matrix.rref();piv=[];col=0
    for i in range(rank):
        while not rr[i,col]:col+=1
        piv.append(col);col+=1
    free=[j for j in range(matrix.ncols()) if j not in set(piv)]
    result=fmpq_mat(matrix.ncols(),len(free))
    for k,j in enumerate(free):
        result[j,k]=1
        for i,p in enumerate(piv):result[p,k]=-rr[i,j]
        divisor=math.lcm(*(rational(result[i,k]).denominator for i in range(result.nrows())))
        values=[int(result[i,k]*divisor) for i in range(result.nrows())]
        common=math.gcd(*values)
        if next(v for v in values if v)<0:common=-common
        for i,v in enumerate(values):result[i,k]=v//common
    assert zero(matrix*result)
    return result,rank

def independent_gram_indices(matrix,prime=1000003):
    n=matrix.nrows();a=nmod_mat(n,n,prime)
    for i in range(n):
        for j in range(n):
            value=rational(matrix[i,j]);assert value.denominator%prime
            a[i,j]=value.numerator*pow(value.denominator,-1,prime)
    rr,rank=a.rref();piv=[];col=0
    for i in range(rank):
        while not rr[i,col]:col+=1
        piv.append(col);col+=1
    sub=nmod_mat([[a[i,j] for j in piv] for i in piv],prime)
    determinant=int(sub.det());assert determinant
    return piv,determinant

def select(matrix,rows,cols):return fmpq_mat([[matrix[i,j] for j in cols] for i in rows])

def matrix_file(path,compute):
    if path.exists():return import_matrix(read(path)['matrix'])
    matrix,profiles=compute()
    write(path,{'matrix':export_matrix(matrix),'profiles':profiles})
    return matrix

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--length',type=int,required=True)
    ap.add_argument('--metric-only',action='store_true');args=ap.parse_args()
    length=args.length;started=time.monotonic()
    centre=read(HERE/f'frames_Y{length}_offset0.json.gz')
    source=read(HERE/f'frames_Y{length}_offset-1.json.gz')
    basis=[r for d,rs in sorted(centre['basis'].items(),key=lambda x:int(x[0])) for r in import_rows(rs)]
    incoming=[r for d,rs in sorted(source['images'].items(),key=lambda x:int(x[0])) for r in import_rows(rs)]
    labels=[(int(d),i) for d,rs in sorted(centre['basis'].items(),key=lambda x:int(x[0])) for i in range(len(rs))]
    print('Y',length,'centre',len(basis),'incoming',len(incoming),flush=True)
    matrix=matrix_file(HERE/f'full_gram_Y{length}.json.gz',lambda:gram(basis+incoming,label=f'Y{length} full'))
    full=select(matrix,range(len(basis)),range(len(basis)))
    independent,det=independent_gram_indices(full)
    expected={5:35,6:379,7:605,8:149}[length]
    assert len(independent)==expected,(length,len(independent),expected)
    physical=select(full,independent,independent)
    # Exact relations are determined by all pairings with a complete independent set.
    embedding=physical.solve(select(full,independent,range(len(basis))))
    assert embedding.transpose()*physical*embedding==full
    relation_indices=[j for j in range(len(basis)) if j not in set(independent)]
    relations=fmpq_mat(len(basis),len(relation_indices))
    for k,j in enumerate(relation_indices):
        relations[j,k]=1
        for i,b in enumerate(independent):relations[b,k]=-embedding[i,j]
    assert zero(full*relations)
    assert len(relation_indices)==len(basis)-expected
    record={'length':length,'full_formal_dimension':len(basis),'finite_dimension':expected,
            'independent_basis_indices':independent,'independent_basis_labels':[labels[j] for j in independent],
            'gram_minor_prime':1000003,'gram_minor_determinant':det,
            'null_relation_indices':relation_indices,'null_relations':export_matrix(relations),
            'embedding_of_all_formal_rows':export_matrix(embedding),'physical_gram':export_matrix(physical),
            'coupling_scope':'Tree radial Wick form and leading interacting differential only; no higher-order operator comparison.'}
    write(HERE/f'finite_metric_Y{length}.json.gz',record)
    print('Y',length,'metric finite rank',expected,'exact nulls',len(relation_indices),'seconds',time.monotonic()-started,flush=True)
    if args.metric_only:return
    qrows=[q(basis[j]) for j in independent]
    outgoing=matrix_file(HERE/f'outgoing_gram_Y{length}.json.gz',lambda:gram(qrows,label=f'Y{length} outgoing'))
    incoming_pairing=select(matrix,range(len(basis),len(basis)+len(incoming)),independent)
    constraints=fmpq_mat([[outgoing[i,j] for j in range(expected)] for i in range(expected)]+
                         [[incoming_pairing[i,j] for j in range(expected)] for i in range(len(incoming))])
    harmonics,rank=rational_nullspace(constraints)
    harmonic_metric=harmonics.transpose()*physical*harmonics
    assert zero(outgoing*harmonics) and zero(incoming_pairing*harmonics)
    write(HERE/f'harmonic_Y{length}.json.gz',{'length':length,'constraint_rank':rank,'harmonic_dimension':harmonics.ncols(),
          'finite_basis_indices':independent,'harmonic_coefficients':export_matrix(harmonics),
          'harmonic_gram':export_matrix(harmonic_metric),'incoming_pairing':export_matrix(incoming_pairing),
          'outgoing_gram':export_matrix(outgoing),'passed':True,'seconds':time.monotonic()-started})
    print('Y',length,'harmonics',harmonics.ncols(),'seconds',time.monotonic()-started,flush=True)

if __name__=='__main__':main()
