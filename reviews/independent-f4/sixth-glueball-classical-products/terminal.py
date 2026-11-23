"""Highest-weight scalar trace polynomials at the first target after removing one glueball."""
from collections import Counter, defaultdict
from functools import cache
from itertools import permutations, product
from pathlib import Path
import json
import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
DEFAULT = HERE
SHAPES = ((2,2,2,2,3), (2,2,2,5), (2,2,3,4), (2,2,7), (2,3,3,3),
          (2,3,6), (2,4,5), (3,3,5), (3,4,4), (4,7), (5,6), (11,), (2,9), (3,8), (11,))
# Base 8 stores a scalar label between zero and six. The sentinel records word length.
ADJOINT = 1 << 35


def sign(order):
    return (-1)**sum(i>j for k,i in enumerate(order) for j in order[k+1:])


@cache
def flags(shape):
    columns = [list(permutations(range(h))) for h in shape]
    signs = [[sign(q) for q in col] for col in columns]
    indices = list(product(*(range(len(c)) for c in columns)))
    digits = np.asarray([sum((c[i] for c,i in zip(columns,idx)),()) for idx in indices], dtype=np.int64).T.copy()
    parity = np.prod(np.array([[s[i] for s,i in zip(signs,idx)] for idx in indices],dtype=np.int64), axis=1)
    return digits,parity


def encode(digits):
    value=np.zeros(digits.shape[1],dtype=np.int64)
    for row in digits:
        value = value*8+row
    return value


def necklace(digits, adjoint):
    n=len(digits)
    def rotations(w):
        raw=encode(w)
        minword=raw.copy()
        for _ in range(n-1):
            raw=((raw<<3)%(8**n))+(raw>>(3*n-3))
            minword=np.minimum(minword,raw)
        return minword
    fwd,rev=rotations(digits), rotations(digits[::-1])
    phase=np.sign(rev-fwd) if n%2 else 1
    return 8**n+ADJOINT*int(adjoint)+np.minimum(fwd,rev),phase


def recipe(shape,seed):
    digits,signs=flags(shape)
    rng=np.random.Generator(np.random.PCG64(926731+seed))
    cuts=SHAPES[seed%len(SHAPES)]
    order=np.split(rng.permutation(11),np.cumsum(cuts)[:-1])
    labels=[]; coefficients=signs.copy()
    for slots in order:
        key,phase=necklace(digits[slots], bool(rng.integers(0,2)))
        labels.append(key);coefficients*=phase
    sorted_rows=np.sort(np.asarray(labels).T,axis=1)
    keep=np.flatnonzero(coefficients)
    if not len(keep):
        return {}
    words,indices=np.unique(sorted_rows[keep],axis=0,return_inverse=True)
    values=np.bincount(indices,weights=coefficients[keep]).astype(np.int64)
    return {tuple(map(int,k)):int(v) for k,v in zip(words,values) if v}


class PolynomialFrame:
    def __init__(self,shape,seeds):
        self.height=shape[0]
        monomials={}
        coeff=[];cols=[];ptr=[0]
        for seed in seeds:
            poly=recipe(shape,seed)
            for word,count in poly.items():
                cols.append(monomials.setdefault(word,len(monomials)));coeff.append(count)
            ptr.append(len(cols))
        self.matrix=csr_matrix((np.asarray(coeff,np.int64),cols,ptr),shape=(len(seeds),len(monomials)))
        atoms=sorted({atom for word in monomials for atom in word})
        self.natoms=len(atoms)
        lookup={k:i for i,k in enumerate(atoms)}
        families=defaultdict(list)
        for w,i in monomials.items():
            families[len(w)].append((i,[lookup[k] for k in w]))
        self.factors=[(np.array([i for i,_ in f]),np.array([idx for _,idx in f])) for f in families.values()]
        self.halves={}
        for ri in (0,1):
            needed={0:{1}};groups=defaultdict(list)
            for i,atom in enumerate(atoms):
                rep,encoded=divmod(atom,ADJOINT)
                if rep!=ri:continue
                n=(encoded.bit_length()-1)//3
                half=n//2
                left=encoded//(8**(n-half));right=8**(n-half)+encoded%(8**(n-half))
                needed.setdefault(half,set()).add(left);needed.setdefault(n-half,set()).add(right)
                groups[half,n-half].append((i,left,right))
            for n in range(max(needed),0,-1):
                needed.setdefault(n-1,set()).update(x//8 for x in needed[n])
            ix={n:{k:i for i,k in enumerate(sorted(vs))} for n,vs in needed.items()}
            levels=[(np.array([ix[n-1][k//8] for k in sorted(needed[n])]),
                     np.array([k%8 for k in sorted(needed[n])])) for n in range(1,max(needed)+1)]
            specs=[(n,m,np.array([i for i,_,_ in row]),np.array([ix[n][j] for _,j,_ in row]),
                    np.array([ix[m][k] for _,_,k in row])) for (n,m),row in groups.items()]
            self.halves[ri]=levels,specs

    def evaluate(self,frame,seed,p):
        rng=np.random.Generator(np.random.PCG64(89223+seed))
        c=rng.integers(-3,4,size=(self.height,52))
        values=np.zeros(self.natoms,np.int64)
        for ri,(levels,specs) in self.halves.items():
            rep=frame[ri];dim=rep.shape[1]
            mats=((c @ rep.reshape(52,-1)).reshape(-1,dim,dim)%p).astype(np.float64)
            assert dim**2 * (p-1)**2 < 2**53
            power=[np.eye(dim)[None]]
            for prefix,last in levels:
                prev=power[-1]
                grp=np.empty((len(last),dim,dim),np.float64)
                for j in range(0,len(last),96):
                    s=slice(j,j+96)
                    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
                        values0=prev[prefix[s]]@mats[last[s]]
                    assert np.isfinite(values0).all() and values0.min()>=0 and values0.max()<=dim*(p-1)**2
                    grp[s]=values0%p
                power.append(grp)
            for a,b,idx,ls,rs in specs:
                for j in range(0,len(idx),64):
                    s=slice(j,j+64)
                    v=np.einsum('nij,nji->n',power[a][ls[s]],power[b][rs[s]])
                    values[idx[s]]=v.astype(np.int64)%p
        products=np.empty(self.matrix.shape[1],np.int64)
        for where,ids in self.factors:
            v=np.ones(len(where),np.int64)
            for atom in ids.T:v=v*values[atom]%p
            products[where]=v
        return np.asarray(self.matrix@products)%p
