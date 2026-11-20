"""Reduce a filtered matrix over a prime field, cancelling each pivot only once."""
import math
import numpy as np
from flint import nmod_mat
from scipy.linalg.blas import dgemm


def pivots(arr,prime):
    if not arr.size or not np.any(arr):
        return []
    reduced,rank=nmod_mat(arr.tolist(),prime).rref()
    result=[]
    col=0
    for i in range(rank):
        while not reduced[i,col]:
            col+=1
        result.append(col)
        col+=1
    return result


def matmul(a,b,prime):
    assert a.shape[1]==b.shape[0]
    out=np.zeros((a.shape[0],b.shape[1]),dtype=np.int64)
    if not out.size:
        return out
    # Each product and every partial dot sum is an integer below 2^53.
    width=min(1024,(2**53-1)//((prime-1)**2))
    assert width>0
    for k in range(0,a.shape[1],width):
        part=dgemm(1.0,a[:,k:k+width],b[k:k+width,:])
        out=(out+part.astype(np.int64)%prime)%prime
    return out


def filtered_ranks(blocks,source_dimensions,target_dimensions,prime):
    assert len(source_dimensions)==len(target_dimensions)
    grades=len(source_dimensions)
    rows=list(source_dimensions)
    cols=list(target_dimensions)
    matrices={key:np.asarray(value%prime,dtype=np.int64) for key,value in blocks.items()}
    counts={}
    for b in range(grades):
        for a in range(b,-1,-1):
            pivot=matrices.get((a,b))
            if pivot is None or not pivot.size:
                counts[a,b]=0
                continue
            assert pivot.shape==(rows[a],cols[b])
            jj=pivots(pivot,prime)
            rank=counts[a,b]=len(jj)
            if not rank:
                continue
            ii=pivots(pivot[:,jj].T,prime)
            assert len(ii)==rank
            square=nmod_mat(pivot[np.ix_(ii,jj)].tolist(),prime)
            left={i:block[:,jj] for (i,j),block in matrices.items() if j==b and i<a}
            right={j:block[ii,:] for (i,j),block in matrices.items() if i==a and j>b}
            # The two edges of the pivot block must also be transformed.
            left[a]=pivot[:,jj]
            right[b]=pivot[ii,:]
            row_mask=np.ones(rows[a],dtype=bool)
            col_mask=np.ones(cols[b],dtype=bool)
            row_mask[ii]=False
            col_mask[jj]=False
            left[a]=left[a][row_mask]
            right[b]=right[b][:,col_mask]
            for (i,j),block in list(matrices.items()):
                if i==a:
                    block=block[row_mask,:]
                if j==b:
                    block=block[:,col_mask]
                matrices[i,j]=block
            for j,edge in right.items():
                if not edge.size or not np.any(edge):
                    continue
                solved=np.array(square.solve(nmod_mat(edge.tolist(),prime)).tolist(),dtype=np.int64)
                for i,entry in left.items():
                    if not entry.size or not np.any(entry):
                        continue
                    product=matmul(entry,solved,prime)
                    matrices[i,j]=(matrices.get((i,j),0)-product)%prime
            assert not np.any(matrices.pop((a,b)))
            rows[a]-=rank
            cols[b]-=rank
    return [[lo,hi,sum(r for (a,b),r in counts.items() if lo<=a<=b<=hi),counts[lo,hi]]
            for lo in range(grades) for hi in range(lo,grades)]
