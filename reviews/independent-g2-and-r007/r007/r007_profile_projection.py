"""Exact orthogonal projection onto all SU(3) boundaries in one predecessor profile."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
import json
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat

import r007_adjoint_matrix_check as generators
import r007_boundary_projection as wick

jets=wick.jets
PROFILE=((0,0,(0,)),)*5+((0,0,(1,)),(0,2,(1,2)),(2,0,(2,)))


def traces():
    lengths=[p for k in range(1,5) for p in combinations_with_replacement(range(2,9),k) if sum(p)==8]
    words=set()
    for sequence in set(permutations(PROFILE)):
        for sizes in lengths:
            off=0
            factors=[]
            for size in sizes:
                factors.append(sequence[off:off+size])
                off+=size
            sign,w=jets.canonical(tuple(factors))
            assert sign==1
            words.add(w)
    return sorted(words)


def coordinates(words):
    # Put the repeated even matrix on a diagonal. The remaining three slots
    # are multilinear, including the odd slot after stripping its sole marker.
    basis=[t for t,_,_ in generators.trace_basis()]
    roots=[tuple(int(h==i)-int(h==j) for h in range(3))
           for i in range(3) for j in range(3) if i!=j]+[(0,0,0)]*2
    triples=[ijk for ijk in product(range(8),repeat=3)
             if all(sum(roots[i][h] for i in ijk)==0 for h in range(3))]
    probes=[(a,b) for a in range(6) for b in triples]
    letters=sorted(set(PROFILE))
    mats={letters[0]:np.array([np.diag([1,t,-1-t]) for t,_ in probes],dtype=np.int64)}
    for k,letter in enumerate(letters[1:]):
        mats[letter]=np.array([basis[s[k]] for _,s in probes],dtype=np.int64)
    @lru_cache(None)
    def trace(word):
        z=mats[word[0]]
        for letter in word[1:]:
            z=z@mats[letter]
        return z.trace(axis1=1,axis2=2)
    cols=[]
    for word in words:
        v=np.ones(len(probes),dtype=np.int64)
        for factor in word:
            v=v*trace(factor)
        cols.append(v)
    # Each nonzero entry of a singleton matrix has magnitude at most 6
    # for the repeated matrix and 2 for any other matrix. Every trace
    # contains at most three diagonal terms at each multiplication of these
    # diagonal or one-entry matrices.
    assert 3**4*6**5*2**3<2**63
    return fmpz_mat(np.array(cols).T.tolist()),len(triples)


def pivots(r,rank):
    return [next(j for j in range(r.ncols()) if r[i,j]) for i in range(rank)]


def coefficients(expr, words):
    return fmpq_mat([[str(expr.get(w,0))] for w in words])


def profile_of(w):
    return tuple(sorted(l for tr in w for l in tr))


def part(state):
    return {w:c for w,c in state.items() if profile_of(w)==PROFILE}


def run():
    target=wick.radial.target()
    k=part(wick.adj.adjoint(target,3))
    words=traces()
    evals,triples=coordinates(words)
    evaluation=fmpq_mat(evals)
    red,dim=evaluation.rref()
    cols=pivots(red,dim)
    small=[words[i] for i in cols]
    @lru_cache(None)
    def adjoint_word(w):
        return tuple(part(wick.adj.adjoint({w:1},3)).items())
    columns=[]
    for w in small:
        q=jets.act({w:1},('Q',0,0))
        v=Counter()
        for term,c in q.items():
            v.update({u:c*d for u,d in adjoint_word(term)})
        columns.append([str(v.get(u,0)) for u in words])
    mat=evaluation*fmpq_mat(columns).transpose()
    rhs=evaluation*coefficients(k,words)
    augmented=fmpq_mat([list(r)+[rhs[i,0]] for i,r in enumerate(mat.tolist())])
    reduced,rank=augmented.rref()
    chosen=pivots(reduced,rank)
    assert all(i<len(small) for i in chosen)
    u={small[j]:Fraction(str(reduced[i,len(small)])) for i,j in enumerate(chosen) if reduced[i,len(small)]}
    boundary=jets.act(u,('Q',0,0))
    assert not jets.act(boundary,('Q',0,0)) and not jets.rho(boundary)
    remainder=Counter(k)
    remainder.subtract(part(wick.adj.adjoint(boundary,3)))
    assert not any((evaluation*coefficients(remainder,words)).entries())
    norm=wick.inner(u,k,3)
    full=Fraction(5171748864)
    single=Fraction(14131200)**2/293253120
    assert norm>=single and Fraction(7987200)<=full-norm<full
    assert wick.inner(boundary,target,3)==norm
    return dict(predecessor_profile=PROFILE,formal_trace_words=len(words),
                torus_neutral_triples=triples,diagonal_interpolation_values=list(range(6)),
                interpolation_coordinates=evals.nrows(),finite_predecessor_dimension=dim,
                boundary_dimension=rank,primitive_terms=len(u),boundary_trace_terms=len(boundary),
                primitive=[[w,str(c)] for w,c in sorted(u.items())],
                projected_adjoint_residual_nonzero_formal_terms=sum(bool(c) for c in remainder.values()),
                projected_adjoint_residual_nonzero_SU3_coordinates=0,
                full_norm=str(full),subtracted_boundary_squared_norm=str(norm),
                squared_norm_after_profile_projection=str(full-norm),
                classical_Cartan_lower_bound='7987200')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    if not __debug__:
        raise RuntimeError('run with assertions enabled')
    record=run()
    args.out.write_text(json.dumps(record,indent=2)+'\n')
    print({k:v for k,v in record.items() if k not in ('primitive','predecessor_profile')},flush=True)
