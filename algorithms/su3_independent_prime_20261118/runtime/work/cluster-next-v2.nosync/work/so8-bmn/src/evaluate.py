"""Exact modular evaluation of both SO(8) invariant parities."""
import numpy as np
from math import comb
from fortuity.hwcoh import evalg
from invariants import ferm

_matrix_size=evalg.matrix_size
_algebra_dim=evalg.algebra_dim
evalg.matrix_size=lambda group,N:8 if group=='SO_odd' and N==8 else _matrix_size(group,N)
evalg.algebra_dim=lambda group,N:28 if group=='SO_odd' and N==8 else _algebra_dim(group,N)


class Frame(evalg.Frame):
    def __init__(self,degree,npts,seed,prime,k=None):
        # SO_odd's antisymmetric generator is reused at size 8. This patch does
        # not reuse any SO_odd dimensions or Cartan restriction.
        super().__init__(8,npts,degree if k is None else k,[degree],seed=seed,p=prime,group='SO_odd')
        self.matrices={};self.scalars={};self.pfaff_cache={};self.pair_cache={}

    def matrix(self,w):
        if w not in self.matrices:
            if len(w)==1:
                d,v=self.L[w[0]]
            else:
                v,d=self.matrix(w[:-1]);v,d=self.mat_mul_letter(v,d,w[-1])
            self.matrices[w]=(v,d)
        return self.matrices[w]

    def trace(self,w):
        if w not in self.scalars:
            v,d=self.matrix(w)
            self.scalars[w]=(np.trace(v,axis1=-2,axis2=-1)%self.p,d)
        return self.scalars[w]

    def color_wedge(self,A,dA,rA,B,dB,rB):
        gi,gj,gs,gst=evalg.table(self.k,dA,dB)
        ci,cj,cs,cst=evalg.table(8,rA,rB)
        R=(A[:,gi[:,None],ci[None,:]]*B[:,gj[:,None],cj[None,:]])%self.p
        R*=gs[None,:,None]
        R*=cs[None,None,:]
        R=np.add.reduceat(R,cst,axis=2)%self.p
        R=np.add.reduceat(R,gst,axis=1)%self.p
        return R

    def epsilon(self,paths):
        # Exterior multiplication of four 2-forms is the epsilon contraction.
        if paths in self.pfaff_cache:return self.pfaff_cache[paths]
        A,d=self.path_pair(paths[:2]);B,e=self.path_pair(paths[2:])
        value=(self.color_wedge(A,d,4,B,e,4)[:,:,0],d+e)
        self.pfaff_cache[paths]=value
        return value

    def path_pair(self,paths):
        if paths not in self.pair_cache:
            rows,cols=np.triu_indices(8,1)
            M,d=self.matrix(paths[0]);N,e=self.matrix(paths[1])
            A=(M[:,:,rows,cols]-M[:,:,cols,rows])%self.p
            B=(N[:,:,rows,cols]-N[:,:,cols,rows])%self.p
            self.pair_cache[paths]=(self.color_wedge(A,d,2,B,e,2),d+e)
        return self.pair_cache[paths]

    def invariant(self,inv):
        paths,traces=inv
        A,d=self.epsilon(paths) if paths else (np.ones((self.P,1),dtype=np.int64),0)
        for w in traces:
            B,e=self.trace(w);A,d=self.scal_mul(A,d,B,e)
        return A.reshape(-1)

    def values(self,words):
        return np.array([self.invariant(w) for w in words],dtype=np.int64).reshape(len(words),self.ncols)
