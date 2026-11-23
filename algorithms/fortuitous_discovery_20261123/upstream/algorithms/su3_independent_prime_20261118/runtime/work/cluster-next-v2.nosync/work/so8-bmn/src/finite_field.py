"""Integer-only FLINT rank. No floating-point tolerance or BLAS kernel."""
import sys,pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'deps'))
from flint import nmod_mat
import numpy as np


def rank_mod(matrix,p):
    matrix=np.asarray(matrix,dtype=np.int64)
    if not matrix.size:return 0
    # Construct row by row to avoid an additional full Python-int matrix.
    m,n=matrix.shape
    result=nmod_mat(m,n,p)
    for i in range(m):
        for k,v in enumerate(matrix[i].tolist()):result[i,k]=v
    return int(result.rank())
