"""The frozen staged elimination with bounded Python-to-FLINT conversion.

No full Python list of matrix entries is formed. FLINT reduces its matrix
in place; output is copied one row at a time to an int32 array.
"""
import sp_environment
import gc
import json
import time

import numpy as np
from flint import nmod_mat
from exact_blas_rows import matmul_mod
from staged_sp2_first_correction import compact


def rref(values,prime):
    assert values.dtype==np.int32 and values.ndim==2
    assert not values.size or (values.min()>=0 and values.max()<prime)
    height,width=values.shape
    matrix=nmod_mat(height,width,prime)
    for i,row in enumerate(values):
        for j,value in enumerate(row):matrix[i,j]=int(value)
    reduced,rank=matrix.rref(True)
    rank=int(rank)
    rows=np.empty((rank,width),dtype=np.int32)
    for i in range(rank):
        rows[i]=np.fromiter((int(reduced[i,j]) for j in range(width)),dtype=np.int32,count=width)
    return rows,rank


def reduce_target(values,columns,previous,variables,prime):
    start=time.monotonic()
    assert values.dtype==np.int32 and columns[-1]==variables
    values,columns=compact(values,columns,variables)
    for number,stage in enumerate(previous):
        lookup={int(c):i for i,c in enumerate(columns)}
        common=[(j,lookup[p]) for j,p in enumerate(stage['pivots']) if p in lookup]
        if not common:continue
        coefficients=values[:,[j for i,j in common]]
        active=np.any(coefficients!=0,axis=0)
        if not np.any(active):continue
        indices=np.flatnonzero(active)
        coefficients=coefficients[:,indices]
        rows=stage['rows'][[common[i][0] for i in indices]]
        union=np.union1d(columns,stage['columns'])
        expanded=np.zeros((len(values),len(union)),dtype=np.int32)
        expanded[:,np.searchsorted(union,columns)]=values
        target_columns=np.searchsorted(union,stage['columns'])
        for begin in range(0,len(values),256):
            end=min(begin+256,len(values))
            product=matmul_mod(coefficients[begin:end],rows,prime)
            expanded[begin:end,target_columns]=(expanded[begin:end,target_columns]-product)%prime
        del values,coefficients,rows
        values,columns=compact(expanded,union,variables)
        del expanded
        assert not set(map(int,columns)).intersection(stage['pivots'])
        print(json.dumps(dict(stage='previous_pivots_eliminated_bounded_flint',previous_stage=number,
                              equations=len(values),columns=len(columns),seconds=time.monotonic()-start)),flush=True)
    values=values[np.any(values!=0,axis=1)]
    if not len(values):
        return dict(rows=np.zeros((0,len(columns)),dtype=np.int32),columns=columns,pivots=[],
                    consistent=True,rank=0,seconds=time.monotonic()-start)
    rows,rank=rref(values,prime)
    del values
    gc.collect()
    local_pivots=[int(np.flatnonzero(row)[0]) for row in rows]
    assert np.array_equal(rows[:,local_pivots],np.eye(rank,dtype=np.int32))
    pivots=[int(columns[i]) for i in local_pivots]
    assert not set(pivots).intersection(p for stage in previous for p in stage['pivots'])
    return dict(rows=rows,columns=columns,pivots=pivots,
                consistent=not pivots or pivots[-1]!=variables,rank=rank,seconds=time.monotonic()-start)
