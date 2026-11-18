"""Modular row extension with provably exact binary64 dot products.

Inputs are integer residues in [0,p). A block has k*(p-1)^2 < 2^53;
every scalar product and every nonnegative partial sum is then an exactly
representable binary64 integer. Reduction is performed in int64.
"""
import sp_environment
import numpy as np
from flint import nmod_mat
from sp_projected_basis import independent_rows


def exact_mod_product(left,right,prime):
    assert left.ndim==right.ndim==2 and left.shape[1]==right.shape[0]
    assert not left.size or (int(left.min())>=0 and int(left.max())<prime)
    assert not right.size or (int(right.min())>=0 and int(right.max())<prime)
    maximum=((1<<53)-1)//((prime-1)**2)
    assert maximum>0
    result=np.zeros((left.shape[0],right.shape[1]),dtype=np.int64)
    for start in range(0,left.shape[1],maximum):
        stop=min(start+maximum,left.shape[1])
        assert (stop-start)*(prime-1)**2 < 2**53
        product=np.asarray(left[:,start:stop],dtype=np.float64)@np.asarray(right[start:stop],dtype=np.float64)
        assert np.isfinite(product).all()
        assert not product.size or (float(product.min())>=0 and float(product.max())<2**53)
        result=(result+product.astype(np.int64))%prime
    return result


def extend_rows(previous,pivots,values,prime):
    free=np.ones(values.shape[1],dtype=bool);free[pivots]=False;free=np.flatnonzero(free)
    reduced=values[:,free].copy()
    if pivots:reduced=(reduced-exact_mod_product(values[:,pivots],previous[:,free],prime))%prime
    ids=independent_rows(reduced,prime)
    if not ids:return previous,pivots,ids
    matrix,rank=nmod_mat(reduced[ids].tolist(),prime).rref();assert int(rank)==len(ids)
    new_free=np.asarray(matrix.tolist(),dtype=np.int64)
    local_pivots=[int(np.flatnonzero(row)[0]) for row in new_free]
    new_pivots=[int(free[j]) for j in local_pivots]
    if len(previous):
        old_free=(previous[:,free]-exact_mod_product(previous[:,new_pivots],new_free,prime))%prime
        previous=previous.copy();previous[:,free]=old_free
    new=np.zeros((len(ids),values.shape[1]),dtype=np.int64);new[:,free]=new_free
    result=np.concatenate((previous,new),axis=0);all_pivots=pivots+new_pivots
    assert np.array_equal(result[:,all_pivots],np.eye(len(result),dtype=np.int64))
    return result,all_pivots,ids
