"""Polarized tree variations with two open words and a colour projector."""
from functools import cache
from math import gcd
import numpy as np
from flag_rows import flags,factor,PAIRS


@cache
def token(slot):
    d=flags()[0]
    if len(slot)==1:return d[slot[0]],1
    a,b=slot;i,j=d[a],d[b]
    low=np.minimum(i,j);high=np.maximum(i,j)
    number=np.zeros(i.shape,dtype=np.int64)
    for k,(u,v) in enumerate(PAIRS):number+=(4+k)*((low==u)&(high==v))
    return number,np.sign(j-i)


def opened(slots,transpose_sign):
    sign=np.ones(flags()[1].shape,dtype=np.int32)
    left=np.ones(sign.shape,dtype=np.int64)
    right=np.full(sign.shape,16**len(slots),dtype=np.int64)
    for j,t in enumerate(slots):
        digits,phase=token(t)
        left=16*left+digits
        right+=16**j*digits
        sign*=phase
    if (-1)**len(slots)*transpose_sign==-1:
        sign*=np.sign(right-left).astype(np.int32)
    return np.minimum(left,right),sign


def choices(y):
    return [(a,b,y-a-b) for a in range(2,7) for b in range(2,a+1)
            if y-a-b in (0,2,3,4,5)]


def polynomial(mode,seed,casimir):
    rng=np.random.default_rng(167819+seed)
    kind='pair' if mode=='out' else ('two','two','three')[seed%3]
    y=10 if mode=='out' else 9
    shapes=choices(y)
    if seed>=10000:
        shapes=[x for x in shapes if x[0]>=4 and x[1]>=4]
    shape=shapes[(seed if mode=='out' else seed//3)%len(shapes)]
    order=tuple(int(i) for i in rng.permutation(11))
    if kind=='three':
        tokens=[order[:3]]+[(i,) for i in order[3:]]
    elif kind=='two':
        tokens=[order[:2],order[2:4]]+[(i,) for i in order[4:]]
    else:tokens=[order[:2]]+[(i,) for i in order[2:]]
    rng.shuffle(tokens)
    terms=[]
    if kind=='three':
        i=next(k for k,t in enumerate(tokens) if len(t)==3)
        a,b,c=tokens[i]
        for x,y0,z in ((a,b,c),(b,c,a),(c,a,b)):
            terms.extend(((1,i,[(x,),(y0,z)]),(-1,i,[(y0,z),(x,)])))
    else:
        for number,i in enumerate(k for k,t in enumerate(tokens) if len(t)==2):
            a,b=tokens[i]
            for sign,pair in ((1,[(a,),(b,)]),(-1,[(b,),(a,)])):
                terms.append((sign*(-1)**number,i,pair))
    keys=[];values=[]
    transpose_sign=1 if casimir==12 else -1
    for phase,i,pair in terms:
        past=0;begin=0;parts=[]
        letters=tokens[:i]+pair+tokens[i+1:]
        sign=flags()[1].copy()*phase
        for h,size in enumerate(shape):
            if not size:
                parts.append(np.zeros(sign.shape,dtype=np.int64));continue
            n=size+(past<=i<past+size)
            section=tuple(letters[begin:begin+n])
            ids,sg=(opened(section,transpose_sign) if h<2 else factor(section,0))
            sign*=sg;parts.append(ids)
            past+=size;begin+=n
        assert begin==len(letters)
        active=sign!=0
        if np.any(active):
            a,b,c=parts
            keys.append(np.stack((np.minimum(a,b)[active],np.maximum(a,b)[active],c[active]),axis=1))
            values.append(sign[active])
    factor.cache_clear()
    if not keys:return {}
    rows,ids=np.unique(np.concatenate(keys),axis=0,return_inverse=True)
    co=np.zeros(len(rows),dtype=np.int64)
    np.add.at(co,ids,np.concatenate(values))
    nz=np.flatnonzero(co)
    d=gcd(*(int(co[i]) for i in nz)) if len(nz) else 1
    return {tuple(map(int,rows[i])):int(co[i])//d for i in nz}
