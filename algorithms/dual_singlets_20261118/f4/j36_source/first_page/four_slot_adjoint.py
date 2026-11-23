"""Contract the full adjoint-matrix colours of the four-, four-, three-, three-slot flag."""

import argparse
from functools import cache
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path

from flint import nmod_mat
import numpy as np
from scipy.sparse import csr_matrix, kron

from colour import Colour
from dotted_filtration import exterior
from flag_boundaries import adjoint
from four_slot_obstruction import PAIRINGS, as_array, mod, trace_pair


CHANNELS = (1,52,324,1053,1274)
LABELS = ((0,0,0,0),(1,0,0,0),(0,0,0,2),(2,0,0,0),(0,1,0,0))


class Adjoint:
    def __init__(self,p):
        c = Colour()
        self.p = p
        self.eye = np.eye(52, dtype=np.int64)
        self.form = np.array([[c.metric.get((i,j),0) for j in range(52)] for i in range(52)],dtype=np.int64)
        self.gens = adjoint()
        self.dual = np.zeros_like(self.form)
        quadratic = np.zeros_like(self.form)
        factors = [csr_matrix(v) for v in self.gens]
        casimir = csr_matrix((52**2,52**2),dtype=np.int64)
        for (i,j),q in c.inverse.items():
            residue = int(q.numerator)*pow(int(q.denominator),-1,p)%p
            self.dual[i,j] = residue
            quadratic += residue * (self.gens[i]@self.gens[j])
            casimir += residue * kron(factors[i],factors[j].T,format='csr')
        assert not np.any((quadratic-18*self.eye)%p)
        self.cross = casimir
        self.cross.data %= p
        self.cross.eliminate_zeros()
        self.trace_generators = csr_matrix(self.gens.transpose(0,2,1).reshape(52,-1))
        self.flat_generators = csr_matrix(self.gens.reshape(52,-1).T)
        self.trace_inverse = self.dual * pow(18,-1,p) % p
        assert not np.any((self.trace_generators @ self.flat_generators).toarray()-18*self.form)

    def casimir(self,x):
        return (36*x-2*(self.cross @ x.ravel()).reshape(52,52)) % self.p

    def split(self,symmetric,skew):
        p = self.p
        scalar = int(symmetric.trace()) * pow(52,-1,p) % p
        trivial = scalar*self.eye
        symmetric = (symmetric-trivial)%p
        small = ((40*symmetric-self.casimir(symmetric))*pow(14,-1,p))%p
        adj = self.trace_inverse @ (self.trace_generators @ skew.ravel() % p) % p
        adj = (self.flat_generators @ adj % p).reshape(52,52)
        return np.stack((trivial,adj,small,(symmetric-small)%p,(skew-adj)%p))

    def columns(self,point):
        p=self.p
        x=np.random.Generator(np.random.PCG64(144336+point)).integers(-3,4,(4,52))
        matrices=(x @ self.gens.reshape(52,-1)).reshape(4,52,52) % p
        brackets={(i,j):matrices[i]@x[j]%p for i,j in combinations(range(4),2)}
        # Endomorphisms obtained by lowering the second colour with the native metric.
        def tensor(a,b):
            return np.outer(a,b) @ self.form
        marked=np.zeros((52,52),dtype=np.int64)
        for perm,s in zip(((0,1,2,3),(0,2,1,3),(0,3,1,2)),(1,-1,1)):
            i,j,k,l=perm
            marked+=s*(tensor(brackets[i,j],brackets[k,l])+tensor(brackets[k,l],brackets[i,j]))
        all_marks=self.split(marked%p,np.zeros_like(marked))
        tail=tensor(x[0],brackets[1,2])-tensor(x[1],brackets[0,2])+tensor(x[2],brackets[0,1])
        reverse=tensor(brackets[1,2],x[0])-tensor(brackets[0,2],x[1])+tensor(brackets[0,1],x[2])
        three=self.split((tail+reverse)%p,(tail-reverse)%p)
        four=np.zeros_like(marked)
        for i in range(4):
            j,k,l=(z for z in range(4) if z!=i)
            scalar=x[j]@self.form@brackets[k,l]%p
            four+=(-1)**i*scalar*matrices[i]
        assert not np.any(all_marks[[0,1,4]]) and not np.any(three[1])
        answer={}
        for k in (2,3): answer['M'+str(CHANNELS[k])]=(all_marks[k],1)
        for k in (0,2,3,4): answer['T'+str(CHANNELS[k])]=(three[k],-1 if k==4 else 1)
        answer['F52']=(four%p,-1)
        if point==0:
            for key,(v,parity) in answer.items():
                assert np.any(v) and not np.any((v.T@self.form-parity*self.form@v)%p)
                if key.endswith('324') or key.endswith('1053'):
                    assert not np.any((self.casimir(v)-(26 if key.endswith('324') else 40)*v)%p)
        return answer


def layout():
    four=('M324','M1053','F52')
    three=('T1','T324','T1053','T1274')
    # At least one four-slot column consists of two marked brackets.
    for a,b in combinations_with_replacement(four,2):
        if a==b=='F52': continue
        for c,d in combinations_with_replacement(three,2):
            for i,j,k,l in PAIRINGS:
                cells=(a,b,c,d)
                yield (cells[i],cells[j]),(cells[k],cells[l])


def evaluate(engine,point,order):
    tensors=engine.columns(point)
    p=engine.p
    @cache
    def join(a,b):
        A,s=tensors[a]
        B,t=tensors[b]
        ab=A@B%p
        ba=B@A%p
        sym=(ab+s*t*ba)*pow(2,-1,p)%p
        return engine.split(sym,(ab-sym)%p)
    names=sorted({pair for test in order for pair in test})
    blocks=np.stack([join(*pair) for pair in names])
    left=blocks.reshape(len(names),5,-1).transpose(1,0,2).astype(np.float64)
    right=blocks.swapaxes(-1,-2).reshape(len(names),5,-1).transpose(1,2,0).astype(np.float64)
    # The complete integer dot products have no rounding at the tested primes.
    assert 52**2*(p-1)**2 < 2**53
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        gram=left@right
    assert np.isfinite(gram).all() and gram.min()>=0 and gram.max()<=52**2*(p-1)**2
    gram=gram.astype(np.int64)%p
    idx={name:i for i,name in enumerate(names)}
    frame=np.array([gram[:,idx[a],idx[b]] for a,b in order]).ravel()
    if point==0:
        for i in range(0,len(order),11):
            a,b=order[i]
            na,nb=idx[a],idx[b]
            assert all(frame[5*i+k]==trace_pair(blocks[na,k],blocks[nb,k])%p for k in range(5))
            # Resolve the original trace pairing in all five inequivalent colours.
            u,v=(tensors[j][0] for j in a)
            x,y=(tensors[j][0] for j in b)
            assert sum(frame[5*i:5*i+5])%p==trace_pair(u@v%p,x@y%p)%p
    return frame


def run(p, start, stop):
    counts=[[exterior(k).get(w,0) for w in LABELS] for k in (3,4)]
    assert counts==[[1,0,1,1,1],[0,1,1,1,0]]
    engine=Adjoint(p)
    order=list(layout())
    assert len(order)==150
    values=[evaluate(engine,i,order) for i in range(start,stop)]
    return np.array(values,dtype=np.uint16)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime',type=int,choices=(1009,1013),required=True)
    parser.add_argument('--start',type=int,default=0)
    parser.add_argument('--stop',type=int,default=703)
    parser.add_argument('--out',required=True,type=Path)
    a=parser.parse_args()
    if not __debug__ or not 0<=a.start<a.stop<=703:
        parser.error('use assertions and a nonempty subset of the 703 nodes')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    np.save(a.out,run(a.prime,a.start,a.stop),allow_pickle=False)
