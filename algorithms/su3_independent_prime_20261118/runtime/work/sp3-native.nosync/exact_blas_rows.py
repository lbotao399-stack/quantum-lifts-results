"""Finite-field row extension using provably exact binary64 dot products.

All inputs are integer residues. Each product and each partial sum is an
integer below 2**53; binary64 multiplication/addition therefore introduces
no rounding. Larger inner dimensions are split into bounded blocks.
"""
import sp_environment
import numpy as np
from flint import nmod_mat
from scipy.linalg.blas import dgemm


def matmul_mod(left,right,prime):
    assert left.shape[1]==right.shape[0]
    assert left.size==0 or (left.min()>=0 and left.max()<prime)
    assert right.size==0 or (right.min()>=0 and right.max()<prime)
    block=(2**53-1)//((prime-1)**2)
    assert block>=1
    result=np.zeros((left.shape[0],right.shape[1]),dtype=np.int64)
    for start in range(0,left.shape[1],block):
        stop=min(start+block,left.shape[1])
        value=dgemm(1.0,left[:,start:stop].astype(np.float64),right[start:stop].astype(np.float64))
        assert (stop-start)*(prime-1)**2<2**53
        assert np.all(np.isfinite(value)) and np.all(value>=0) and np.all(value<2**53)
        result=(result+value.astype(np.int64)%prime)%prime
    return result


def independent_rows(values,prime):
    if not values.shape[0]:return []
    reduced,rank=nmod_mat(values.T.tolist(),prime).rref()
    return [next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]


def extend_rows(previous,pivots,values,prime):
    reduced=(values-matmul_mod(values[:,pivots],previous,prime))%prime if pivots else values.copy()
    ids=independent_rows(reduced,prime)
    if not ids:return previous,pivots,ids
    matrix,rank=nmod_mat(reduced[ids].tolist(),prime).rref();assert int(rank)==len(ids)
    new=np.array(matrix.tolist(),dtype=np.int64);new_pivots=[int(np.flatnonzero(row)[0]) for row in new]
    assert not set(pivots).intersection(new_pivots)
    if len(previous):previous=(previous-matmul_mod(previous[:,new_pivots],new,prime))%prime
    result=np.concatenate((previous,new),axis=0);all_pivots=pivots+new_pivots
    assert np.array_equal(result[:,all_pivots],np.eye(len(result),dtype=np.int64))
    return result,all_pivots,ids
