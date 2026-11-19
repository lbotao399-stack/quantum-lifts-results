"""A nonzero minor selects words and monomials; equality is checked in every
integer coefficient after rational solution with the original signed entries.
The sampled ranks are lower bounds only. No full-rank assumption is needed.
"""
from fractions import Fraction
from functools import reduce
from math import gcd,lcm
from pathlib import Path
import gzip, json, hashlib
import numpy as np
from flint import nmod_mat, fmpz_mat
from three_content_physical import PRIME, s

def canonical_digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def pivots(a):
    x,r=nmod_mat((a%PRIME).tolist(),PRIME).rref(); cols=[]; j=0
    for i in range(r):
        while not x[i,j]:j+=1
        cols.append(j);j+=1
    return cols

def product(A,rows,coeff,scale):
    # Compute against all monomials with a certified accumulator bound or Python integers.
    x=A[:-1,:].tocsr(); b=A.getrow(A.shape[0]-1).toarray().ravel(); c=[0]*x.shape[0]
    for i,t in zip(rows,coeff):c[i]=int(t)
    bound=max(abs(int(t)) for t in x.data)*sum(abs(v) for v in c)+max(abs(int(v)) for v in b)*abs(scale)
    if bound<2**63-1:
        d=np.asarray(np.asarray(c,dtype=np.int64)@x).ravel()-int(scale)*b;return np.flatnonzero(d).tolist()
    v=[-int(t)*scale for t in b];
    for i in rows:
        a,z=x.indptr[i:i+2]; m=c[i]
        for k,t in zip(x.indices[a:z],x.data[a:z]):v[k]+=m*int(t)
    return [k for k,c in enumerate(v) if c]

def verify(w,z,A,keys):
    result=s.read(w/'solution.json.gz');assert result['rhs_sha256']==z['rhs_sha256'] and result['contents']==z['contents'] and (result['recipe_canonical_sha256']==canonical_digest(z) if 'recipe_canonical_sha256' in result else result['recipe_sha256']==hashlib.sha256((w/'recipe.json.gz').read_bytes()).hexdigest()); sol=result['solution'];assert len(set(tuple(k) for k,_ in sol))==len(sol);indices={tuple(v):i for i,v in enumerate(z['source'])}; rows=[indices[tuple(v)] for v,k in sol];co=[Fraction(k)*z['den'] for v,k in sol];den=lcm(*(v.denominator for v in co));rem=product(A,rows,[int(v*den) for v in co],den); assert not rem and result['verified_polynomial_coordinates']==A.shape[1];print('all polynomial coefficients verified:',A.shape[1],flush=True)

def solve(w,z,A,keys):
    n=A.shape[0]-1;m=A.shape[1];rng=np.random.default_rng(714293);order=rng.permutation(m).tolist();x=A.tocsr();# exact selected three-colour polynomials
    prefer=list(range(n)); samples=[];width=512;mandatory=[]
    while True:
        sample=list(dict.fromkeys(mandatory+order[:width])); F=A[:,sample].toarray();cols=pivots(F[prefer,:]);rank=len(cols);rows0=pivots(F[prefer,:][:,cols].T);rows=[prefer[j] for j in rows0];I=[sample[i] for i in cols];assert len(rows)==rank
        aug=len(pivots(F));print('sample width, ranks:',len(sample),rank,aug,flush=True);samples.append(dict(width=len(sample),rank=rank,augmented_rank=aug))
        if aug!=rank:
            s.write(w/'restricted_gap.json',dict(samples=samples,source_words=n,monomials=m,scope='Modular sampled rank gap on these selected parents only.'));return
        a=F[rows,:][:,cols].T;b=F[-1,cols];M=fmpz_mat(a.tolist());B=fmpz_mat([[int(v)] for v in b]);print('rational solve',rank,flush=True);v=M.solve(B);coef=[Fraction(str(v[i,0])) for i in range(rank)];den=lcm(*(c.denominator for c in coef));non=product(A,rows,[int(c*den) for c in coef],den);print('nonzero polynomial residual',len(non),'of',m,flush=True)
        if not non:
            recipe=hashlib.sha256((w/'recipe.json.gz').read_bytes()).hexdigest();out=dict(rhs_sha256=z['rhs_sha256'],recipe_sha256=recipe,recipe_canonical_sha256=canonical_digest(z),contents=z['contents'],monomial_key_order_sha256=hashlib.sha256(keys.astype('<i8').tobytes()).hexdigest(),source_words=n,trace_words=len(z['word']),verified_polynomial_coordinates=m,rhs_nonzero_coordinates=int(x[-1].nnz),sample_ranks=samples,selected_minor_modulus=PRIME,selected_minor_determinant=int(nmod_mat((a%PRIME).tolist(),PRIME).det()),pivot_word_indices=rows,pivot_monomial_keys=[int(keys[i]) for i in I],solution=[[z['source'][i],str(c/z['den'])] for i,c in zip(rows,coef) if c],nonzero_polynomial_residual=0,scope=z.get('scope','Three specified complete contents only; both side components of these parents are formally zero.'))
            s.write(w/'solution.json.gz',out);verify(w,z,A,keys);return
        # A surviving polynomial coefficient supplies an additional independent condition
        # over Q. Also add random columns in case an omitted direction vanishes mod PRIME.
        mandatory=list(dict.fromkeys(mandatory+non[:16]));assert width<min(2400,m),'restricted solve unresolved';width=min(2400,m,width+256)
