"""Incremental RREF without repeated conversion of the retained matrix to Python lists."""
import sp_environment
import numpy as np
from flint import nmod_mat
from exact_blas_rows import matmul_mod as product


class RowSpace:
    def __init__(self,columns,prime):
        self.columns,self.prime=columns,prime
        self.rows=np.empty((0,columns),dtype=np.int64)
        self.pivots=[]

    def extend(self,values):
        if not len(values):return []
        values=np.asarray(values,dtype=np.int64)%self.prime
        if self.pivots:
            values=(values-product(values[:,self.pivots],self.rows,self.prime))%self.prime
        matrix=nmod_mat(values.tolist(),self.prime)
        transposed,rank=matrix.transpose().rref()
        if not rank:return []
        selected=[next(j for j in range(transposed.ncols()) if transposed[i,j]) for i in range(rank)]
        reduced,_=matrix.rref()
        newrows=np.asarray([[int(reduced[i,j]) for j in range(self.columns)] for i in range(rank)],dtype=np.int64)
        newpivots=[int(np.flatnonzero(row)[0]) for row in newrows]
        if self.pivots:
            self.rows=(self.rows-product(self.rows[:,newpivots],newrows,self.prime))%self.prime
        self.rows=np.concatenate((self.rows,newrows),axis=0)
        self.pivots.extend(newpivots)
        # Pivot order need not be sorted: rows[:,pivots] is exactly the identity.
        assert len(set(self.pivots))==len(self.pivots)
        return selected

    def dump(self):
        return dict(columns=self.columns,prime=self.prime,pivots=self.pivots,rows=self.rows.astype(np.uint32))

    @classmethod
    def restore(cls,state):
        output=cls(state['columns'],state['prime'])
        output.pivots=list(state['pivots'])
        output.rows=np.asarray(state['rows'],dtype=np.int64).reshape(len(output.pivots),output.columns)
        return output
