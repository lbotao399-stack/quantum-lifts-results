"""Integral colour probes in the derivation algebra of the real Albert algebra."""
from functools import lru_cache
from itertools import combinations
import numpy as np
from flint import nmod_mat

import mixed_incoming as m


def octonions():
    table = {(0,i):(1,i) for i in range(8)}
    table.update({(i,0):(1,i) for i in range(8)})
    for i in range(1,8):
        table[i,i] = -1,0
    for i,j,k,s in ((1,2,3,1),(1,4,5,1),(1,6,7,1),(2,4,6,1),
                    (2,5,7,-1),(3,4,7,-1),(3,5,6,-1)):
        for a,b,c in ((i,j,k),(j,k,i),(k,i,j)):
            table[a,b] = s,c
            table[b,a] = -s,c
    assert len(table) == 64
    return table


@lru_cache(None)
def basis():
    table = octonions()
    pairs = ((0,1),(0,2),(1,2))
    labels = [(i,i,0) for i in range(3)]+[(i,j,u) for i,j in pairs for u in range(8)]
    herm = [[(i,j,u,1)] if i==j else [(i,j,u,1),(j,i,u,1 if u==0 else -1)]
            for i,j,u in labels]
    left = np.zeros((27,27,27), dtype=np.int64)
    for i,a in enumerate(herm):
        for j,b in enumerate(herm):
            matrix = np.zeros((3,3,8),dtype=np.int64)
            for x,y in ((a,b),(b,a)):
                for h,k,u,c in x:
                    for kk,l,v,d in y:
                        if k==kk:
                            sg,w = table[u,v]
                            matrix[h,l,w] += c*d*sg
            assert np.array_equal(matrix[:,:,0],matrix[:,:,0].T)
            assert np.array_equal(matrix[:,:,1:],-matrix[:,:,1:].swapaxes(0,1))
            left[i,:,j] = [matrix[label] for label in labels]
    generators = np.asarray([left[i]@left[j]-left[j]@left[i] for i,j in combinations(range(27),2)])
    flat = generators.reshape(-1,27*27)
    ids = m.pivotcols(nmod_mat(flat.T.tolist(),1000003))
    assert len(ids) == 52
    mats = generators[ids]
    unit = np.array([1,1,1]+[0]*24)
    metric = np.diag([1,1,1]+[2]*24)
    assert not (mats@unit).any()
    assert not (mats.swapaxes(-1,-2)@metric + metric@mats).any()
    for d in mats:
        # [D,L_x]=L_(Dx) checks the derivation identity on all pairs of basis elements.
        assert np.array_equal(d@left-left@d, np.einsum("hi,hkj->ikj",d,left))
    return mats,ids


def frame(letters, points, k, degree, prime, seed):
    mats,_ = basis()
    letters = sorted(letters)
    f = m.Frame(27, points, k, [degree], letters, seed=1, p=prime)
    rng = np.random.default_rng(seed)
    for l in letters:
        parity, template = f.L[l]
        c = rng.integers(0,prime,size=template.shape[:-2]+(52,),dtype=np.int64)
        f.L[l] = parity,(c @ mats.reshape(52,-1)).reshape(template.shape) % prime
    return f
