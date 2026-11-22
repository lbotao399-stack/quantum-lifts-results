"""Exact highest-weight cochains in the two level-50 BMN source columns."""
from itertools import permutations
from math import comb
import json
from pathlib import Path
import numpy as np


def characters(weights):
    d=max(map(sum,weights));shape=tuple(max(w[n] for w in weights)+2-n+1 for n in range(3))+(d+1,2*d+1)
    a=np.zeros(shape,dtype=np.int64);a[0,0,0,0,d]=1
    # Twelve commuting and nine anticommuting colour components. This positive
    # bound covers every intermediate coefficient, including the unrefined one.
    assert 2**9*comb(d+11,11)<2**62
    for mask in range(1,8):
        odd=[int(bool(mask&(1<<k))) for k in range(3)]
        for z in (-1,0,1):
            power=1
            while power<=d:
                move=tuple(power*t for t in (*odd,1,z));dst=[];src=[]
                if any(abs(x)>=n for x,n in zip(move,shape)): break
                for step,n in zip(move,shape):
                    dst.append(slice(max(step,0),n+min(step,0)));src.append(slice(max(-step,0),n-max(step,0)))
                a[tuple(dst)]+=a[tuple(src)].copy()
                if mask.bit_count()%2==0: break
                power*=2
    sing=a[:,:,:,:,d]-a[:,:,:,:,d+1]; out={}
    for w in weights:
        dims=np.zeros(d+1,dtype=np.int64)
        for p in permutations(range(3)):
            r=tuple(w[i]+p[i]-i for i in range(3));inv=sum(v>u for i,v in enumerate(p) for u in p[i+1:]);
            if min(r)>=0:dims+=(-1)**inv*sing[r]
        assert (dims>=0).all();out[','.join(map(str,w))]={str(y):int(n) for y,n in enumerate(dims) if n}
    return out

if __name__=='__main__':
    result=characters([(9,8,8),(9,9,7),(4,4,4),(5,4,4),(6,6,6),(7,7,7)]);print(json.dumps(result,indent=2))
