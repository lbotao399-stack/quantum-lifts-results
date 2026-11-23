"""Paired projected words as classical boundaries in the Taylor diagram (5,4,2)."""
from collections import defaultdict
from pathlib import Path
import json
import numpy as np
from scipy.sparse import csr_matrix

from terminal import flags,encode,necklace,DEFAULT


def open_word(digits):
    n=len(digits);u=encode(digits);v=encode(digits[::-1])
    return (8**n+np.minimum(u,v), np.sign(v-u) if n%2 else 1)


def recipe(shape,seed):
    lengths=((4,7),(5,6),(4,5,2))[seed%3]
    rng=np.random.Generator(np.random.PCG64(412771+seed))
    bits,signs=flags(shape)
    indices=np.split(rng.permutation(11),np.cumsum(lengths)[:-1])
    (a,s),(b,t)=[open_word(bits[x]) for x in indices[:2]]
    c,last=necklace(bits[indices[2]],False) if len(lengths)==3 else (np.ones(len(signs),np.int64),1)
    signs=signs*s*t*last
    keep=np.flatnonzero(signs)
    if not len(keep):return {}
    keys,ix=np.unique(np.stack([np.minimum(a,b),np.maximum(a,b),c],axis=1)[keep],axis=0,return_inverse=True)
    values=np.bincount(ix,weights=signs[keep]).astype(np.int64)
    return {tuple(map(int,k)):int(v) for k,v in zip(keys,values) if v}


class Frame:
    def __init__(self,shape,seeds):
        self.h=shape[0]
        where={};co=[];ix=[];ptr=[0]
        for s in seeds:
            for k,v in recipe(shape,s).items():
                ix.append(where.setdefault(k,len(where)));co.append(v)
            ptr.append(len(co))
        self.matrix=csr_matrix((np.asarray(co,np.int64),ix,ptr),shape=(len(seeds),len(where)))
        keys=list(where)
        opened=sorted({w for a,b,c in keys for w in (a,b)})
        closed=sorted({c for a,b,c in keys})
        byword={k:i for i,k in enumerate(opened)};bytrace={k:i for i,k in enumerate(closed)}
        self.a,self.b,self.c=np.array([(byword[a],byword[b],bytrace[c]) for a,b,c in keys]).T
        self.closed=closed
        groups=defaultdict(list);needed={0:{1}}
        for i,w in enumerate(opened):
            n=(w.bit_length()-1)//3;m=n//2
            a=w//(8**(n-m));b=8**(n-m)+w%(8**(n-m))
            needed.setdefault(m,set()).add(a);needed.setdefault(n-m,set()).add(b)
            groups[m,n-m].append((i,a,b))
        for n in range(max(needed),0,-1):
            needed.setdefault(n-1,set()).update(w//8 for w in needed[n])
        idx={n:{w:i for i,w in enumerate(sorted(g))} for n,g in needed.items()}
        self.levels=[(np.array([idx[n-1][w//8] for w in sorted(needed[n])]),
                      np.array([w%8 for w in sorted(needed[n])])) for n in range(1,max(needed)+1)]
        self.groups=[(m,n,np.array([(i,idx[m][a],idx[n][b]) for i,a,b in raw]).T) for (m,n),raw in groups.items()]
        self.size=len(opened)

    def evaluate(self,rep,jordan,gram,point,p):
        c=np.random.Generator(np.random.PCG64(89223+point)).integers(-3,4,size=(self.h,52))
        x=((c@rep.reshape(52,-1)).reshape(-1,26,26)%p).astype(np.float64)
        powers=[np.eye(26)[None]]
        def multiply(a,b):
            with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
                value=a@b
            assert np.all(np.isfinite(value)) and value.min()>=0 and value.max()<=26*(p-1)**2
            return value%p
        for before,last in self.levels:
            powers.append(multiply(powers[-1][before],x[last]))
        projections=np.empty((self.size,26),np.int64)
        for m,n,(ids,aa,bb) in self.groups:
            for i in range(0,len(aa),90):
                s=slice(i,i+90)
                mat=multiply(powers[m][aa[s]],powers[n][bb[s]])
                projections[ids[s]]=np.asarray(jordan@mat.reshape(-1,676).T).T.astype(np.int64)%p
        cov=(projections@gram)%p
        traces=np.array([1 if w==1 else int(np.einsum('ij,ji',x[(w//8)%8],x[w%8]))%p for w in self.closed])
        values=np.empty(len(self.a),np.int64)
        for i in range(0,len(values),4096):
            s=slice(i,i+4096)
            values[s]=(np.einsum('ij,ij->i',cov[self.a[s]],projections[self.b[s]])%p)*traces[self.c[s]]%p
        return self.matrix@values%p
