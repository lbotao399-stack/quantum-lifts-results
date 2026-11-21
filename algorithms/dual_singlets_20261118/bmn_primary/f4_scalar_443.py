"""Exact boundary insertions and two-slot scalar candidates in the F4 (4,4,3) sector."""
import argparse
from itertools import permutations
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import mixed_incoming as m
from f4_colour import frame
from f4_scalar_nine import act, on_traceless, pairing, tensors


def apply(matrix, tensor, axis, prime):
    assert 26*(prime-1)**2 < 2**63
    return act(matrix,tensor,axis,prime)


def alternated(matrices,tensor,prime):
    value=np.zeros_like(tensor)
    for perm in permutations(range(3)):
        sign=(-1)**sum(i>j for k,i in enumerate(perm) for j in perm[k+1:])
        term=tensor
        for axis,matrix in enumerate(perm):
            term=apply(matrices[matrix],term,axis,prime)
        value+=sign*term
    return value%prime


def compute(work, points, prime, seed):
    work.mkdir(parents=True,exist_ok=True)
    letters=[(0,0,(i,)) for i in range(3)]
    f=frame(letters,points,1,0,prime,seed)
    blocks=m.words_by_shape(11,(0,0,4,4,3))
    terms=[e for _,block in sorted(blocks.items()) for e in m.singlets(block)]
    words,coeff=m.flatten(terms)
    raw=m.values(words,f,m.library(work,27,prime))
    baseline=np.zeros((len(coeff),points),dtype=np.int64)
    for i,row in enumerate(coeff):
        for k,n in row:
            baseline[i]=(baseline[i]+(n%prime)*raw[k])%prime
    data=tensors(prime)
    images=[]
    candidates=[]
    for i in range(points):
        scalars=on_traceless([f.L[l][1][i,0] for l in letters],data,prime)
        comm=(scalars[0]@scalars[1]-scalars[1]@scalars[0])%prime
        left=[data[2]]
        right=[data[3]]
        for _ in range(3):
            left.append(alternated([s.T for s in scalars],left[-1],prime))
            right.append(alternated(scalars,right[-1],prime))
        bound=[]
        probes=[]
        for n in range(4):
            v=right[3-n]
            for axis in range(4):
                bound.append(pairing(left[n],apply(comm,v,axis,prime),prime))
            l_first=[apply(scalars[0].T,left[n],s,prime) for s in range(3)]
            l_second=[apply(scalars[1].T,left[n],s,prime) for s in range(3)]
            r_first=[apply(scalars[0],v,s,prime) for s in range(1,4)]
            r_second=[apply(scalars[1],v,s,prime) for s in range(1,4)]
            for a in range(3):
                for b in range(a+1,4):
                    probes.append((pairing(l_first[a],r_second[b-1],prime)-
                                   pairing(l_second[a],r_first[b-1],prime))%prime)
        # Modulo scalar commutators the alternating triple and the alternating
        # pair commute on the tensor product. Their insertion orders differ by Q0 boundaries.
        bound.extend((probes[6*n+s]-probes[s])%prime for n in range(1,4) for s in range(6))
        images.append(bound)
        candidates.append(probes)
    extended=np.vstack((baseline,np.array(images,dtype=np.int64).T))
    combined=np.vstack((extended,np.array(candidates,dtype=np.int64).T))
    result={"prime":prime,"seed":seed,"points":points,"weight":[4,4,3],"length":11,
            "formal_trace_highest_weights":len(terms),"cochain_dimension":14}
    for name,matrix in (("trace",baseline),("boundaries",extended),("with_scalars",combined)):
        cert=m.minor(matrix,prime)
        assert cert["rank"]<=14
        result[name]={"minor":cert,"values":matrix[np.ix_(cert["rows"],cert["columns"])].tolist()}
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--prime",type=int,choices=(1000003,1000033),default=1000003)
    p.add_argument("--seed",type=int,default=834057)
    p.add_argument("--points",type=int,default=15)
    a=p.parse_args()
    d=compute(a.work,a.points,a.prime,a.seed)
    a.out.write_text(json.dumps(d,separators=(",",":"))+"\n")
    print({k:d[k]["minor"]["rank"] for k in ("trace","boundaries","with_scalars")})
