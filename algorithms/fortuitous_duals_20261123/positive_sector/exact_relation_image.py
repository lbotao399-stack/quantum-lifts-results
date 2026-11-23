#!/usr/bin/env python3
"""Exact image-rank certificate for the complete rational identity space."""
from probe_stable import HERE
from fractions import Fraction
from mixed33_exact import exact,apply
from flint import fmpz_mat
from scipy.sparse import load_npz
import json,gzip,math,time
import numpy as np

def main():
    start=time.monotonic()
    data=json.loads(gzip.decompress((HERE/'rational_relation_basis.json.gz').read_bytes()))
    coefficients=np.zeros((len(data['basis_coefficients']),data['formal_dimension']),dtype=np.int64)
    for i,(row,denominator) in enumerate(zip(data['basis_coefficients'],data['integer_row_denominators'],strict=True)):
        for j,c in row:coefficients[i,j]=int(Fraction(c)*denominator)
    q=load_npz(HERE/'stable_q_columns.npz')
    bound=int(abs(q).sum(axis=1).max())*int(abs(coefficients).max())
    assert bound<2**63
    image=np.asarray(q@coefficients.T).T
    rows=[{int(j):int(row[j]) for j in np.flatnonzero(row)} for row in image]
    independent,pivots,relations=exact(rows)
    assert all(not apply(r,rows) for r in relations)
    assert len(set(max(r) for r in relations))==len(relations)
    assert len(independent)+len(relations)==len(rows)
    matrix=fmpz_mat([[rows[i].get(j,0) for j in pivots] for i in independent])
    determinant=int(matrix.det());assert determinant
    assert len(independent)==116 and len(relations)==65
    certificate={'passed':True,'coefficient_field':'Q','image_rank':len(independent),'nullity':len(relations),'input_relations':len(rows),'formal_output_coordinates':q.shape[0],'integer_minor_determinant':str(determinant),'minor_rows':independent,'minor_columns':pivots,'independent_integer_null_relations':[[[i,str(c)] for i,c in sorted(row.items())] for row in relations],'all65_exact_null_residuals_zero':True,'integer_image_absolute_bound':bound,'seconds':time.monotonic()-start}
    (HERE/'exact_relation_image_result.json').write_text(json.dumps(certificate,indent=2)+'\n')
    print(json.dumps({k:v for k,v in certificate.items() if k not in ['minor_rows','minor_columns','independent_integer_null_relations','integer_minor_determinant']}),flush=True)

if __name__=='__main__':main()
