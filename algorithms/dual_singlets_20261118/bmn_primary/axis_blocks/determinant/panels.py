"""Rank a sparse integral matrix using bounded prime-field BLAS panels."""
from ctypes import CDLL, POINTER, c_float, c_int
from pathlib import Path
import numpy as np

def rank(columns, prime):
    if prime not in (127,131): raise ValueError(prime)
    active_rows=sorted({i for col in columns for i,n in col.items() if n % prime})
    active_cols=[col for col in columns if any(n % prime for n in col.values())]
    if not active_rows or not active_cols:return 0
    ids={w:i for i,w in enumerate(active_rows)}
    a,b=len(ids),len(active_cols)
    transposed=a<b
    mat=np.zeros((b,a) if transposed else (a,b),dtype=np.float32)
    for j,col in enumerate(active_cols):
        for i,n in col.items():
            n%=prime
            if n:
                if transposed: mat[j,ids[i]]=n
                else: mat[ids[i],j]=n
    fn=CDLL(str(Path(__file__).with_name('panel_rank.so').resolve())).panel_rank
    fn.argtypes=(POINTER(c_float),c_int,c_int,c_int)
    fn.restype=c_int
    answer=fn(mat.ctypes.data_as(POINTER(c_float)),*mat.shape,prime)
    assert 0<=answer<=min(mat.shape)
    return answer
