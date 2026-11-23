"""Classical boundary rows from four Albert strands in the (5,4,2) diagram."""
from functools import cache
from itertools import combinations,product
from pathlib import Path
import json
import numpy as np

from terminal import DEFAULT

# The first two blocks alternate five and four scalar labels.
FIVE=((0,0,1,2,3),(0,0,1,1,2),(0,0,0,1,2),(2,2,0,1,3),(1,1,2,2,3),(2,2,2,0,3))
FOUR=((0,0,1,2),(1,1,2,3),(1,1,0,2),(2,2,0,3),(0,0,1,3),(1,1,0,3))


def parity(seq):
    return (-1)**sum(i>j for k,i in enumerate(seq) for j in seq[k+1:])


def act(m,tensor,axis):
    aligned=np.ascontiguousarray(np.moveaxis(tensor,axis,0)).reshape(26,-1)
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        result=m@aligned
    return np.moveaxis(result.reshape((26,)*4),0,axis)


def point(number,rep,low,high,p,recipes):
    coeff=np.random.Generator(np.random.PCG64(89223+number)).integers(-3,4,size=(5,52))
    x=(coeff@rep.reshape(52,-1)).reshape(5,26,26)%p
    @cache
    def alternating(labels):
        if not labels:return np.eye(26,dtype=np.int64)
        return sum((-1)**j*x[i]@alternating(labels[:j]+labels[j+1:]) for j,i in enumerate(labels))%p
    @cache
    def end(n,layout,dual):
        root=low if dual else high
        arms=tuple(sorted(set(layout),key=lambda i:(-layout.count(i),i)))
        out=np.zeros((26,)*4,np.float64)
        labels=tuple(range(n))
        def fill(rest,idx,t,assigned):
            if idx==len(arms):
                out[:] += parity(assigned)*t
                return
            axis=arms[idx];size=layout.count(axis)
            for keep in combinations(rest,size):
                complement=tuple(i for i in rest if i not in keep)
                op=alternating(keep).astype(np.float64)
                y=act(op.T if dual else op,t,axis)
                # Each reduction is exact with dimension 26 and p<=1013.
                assert np.all(np.isfinite(y)) and np.abs(y).max()<2**53
                y=y.astype(np.int64)%p
                fill(complement,idx+1,y.astype(np.float64),assigned+keep)
        fill(labels,0,root,())
        assert np.all(np.isfinite(out)) and np.abs(out).max()<2**53
        return out.astype(np.int64)%p
    @cache
    def joined(n,layout,dual,axis):
        v=end(n,layout,dual).astype(np.float64)
        return tuple(np.ascontiguousarray(act(m.T if dual else m,v,axis)).ravel().astype(np.int64)%p for m in x[:2].astype(np.float64))
    result=[]
    for f,g,i,j in recipes:
        a=joined(5,FIVE[f],True,i);b=joined(4,FOUR[g],False,j)
        # Reduce middle actions before the length-26^4 scalar contraction.
        left=a
        right=b
        result.append(int((left[0]@right[1]-left[1]@right[0])%p))
    for fn in (joined,end,alternating):fn.cache_clear()
    return result
