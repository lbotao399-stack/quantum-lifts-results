"""Exact row extension without multiplying through existing identity columns."""
import sp_environment
import numpy as np
from flint import nmod_mat
from sp_projected_basis import independent_rows

def extend_rows(previous, pivots, values, prime):
    assert max(previous.shape[0], values.shape[0])*(prime-1)**2 < 2**63
    free=np.ones(values.shape[1],dtype=bool);free[pivots]=False
    free=np.flatnonzero(free)
    reduced=values[:,free].copy()
    if pivots:
        reduced=(reduced-values[:,pivots]@previous[:,free])%prime
    ids=independent_rows(reduced,prime)
    if not ids:
        return previous,pivots,ids
    matrix,rank=nmod_mat(reduced[ids].tolist(),prime).rref()
    assert int(rank)==len(ids)
    new_free=np.asarray(matrix.tolist(),dtype=np.int64)
    local_pivots=[int(np.flatnonzero(row)[0]) for row in new_free]
    new_pivots=[int(free[j]) for j in local_pivots]
    if len(previous):
        old_free=(previous[:,free]-previous[:,new_pivots]@new_free)%prime
        previous=previous.copy();previous[:,free]=old_free
    new=np.zeros((len(ids),values.shape[1]),dtype=np.int64);new[:,free]=new_free
    result=np.concatenate((previous,new),axis=0)
    all_pivots=pivots+new_pivots
    assert np.array_equal(result[:,all_pivots],np.eye(len(result),dtype=np.int64))
    return result,all_pivots,ids
