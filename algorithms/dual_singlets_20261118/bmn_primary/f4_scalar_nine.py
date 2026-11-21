"""An exceptional scalar singlet from the invariant cubic of the traceless Albert algebra."""
import argparse
from itertools import permutations
import json
from pathlib import Path
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import mixed_incoming as m
from f4_colour import frame, jordan_multiplication


def tensors(prime):
    metric=np.array([1,1,1]+[2]*24,dtype=np.int64)
    embed=np.zeros((27,26),dtype=np.int64)
    embed[:3,0]=(1,-1,0)
    embed[:3,1]=(1,1,-2)
    embed[3:,2:]=np.eye(24,dtype=np.int64)
    gram=(embed.T*metric)@embed
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    inverse=np.array([pow(int(v),-1,prime) for v in gram.diagonal()],dtype=np.int64)
    project=(inverse[:,None]*embed.T*metric)%prime
    assert np.array_equal(project@embed%prime,np.eye(26,dtype=np.int64))
    left=np.einsum("ui,uxy->ixy",embed,jordan_multiplication())
    cubic=np.einsum("xj,x,ixy,yk->ijk",embed,metric,left,embed,optimize=True)
    assert np.array_equal(cubic,cubic.transpose(1,2,0)) and np.array_equal(cubic,cubic.swapaxes(0,1))
    lower=((cubic.reshape(-1,26)*inverse)@cubic.reshape(-1,26).T%prime).reshape((26,)*4)
    upper=lower.transpose(0,2,1,3).copy()
    for axis in range(4):
        shape=[1]*4
        shape[axis]=26
        upper=upper*inverse.reshape(shape)%prime
    return embed,project,lower,upper


def on_traceless(colour, tensors_data, prime):
    embed,project,_,_=tensors_data
    if 27*max(prime-1,6)**2 >= 2**63:
        raise ValueError("A tensor contraction exceeds the signed integer accumulator")
    return [project@x%prime@embed%prime for x in colour]


def act(matrix, value, axis, prime):
    source=np.moveaxis(value,axis,0)
    applied=np.einsum("ij,jklm->iklm",matrix,source,dtype=np.int64,optimize=False)%prime
    return np.moveaxis(applied,0,axis)


def alternate(matrices, value, prime):
    result=np.zeros_like(value)
    for p in permutations(range(3)):
        sign=(-1)**sum(a>b for j,a in enumerate(p) for b in p[j+1:])
        word=value
        for axis,number in enumerate(p):
            word=act(matrices[number],word,axis,prime)
        result+=sign*word
    return result%prime


def pairing(lower, upper, prime):
    return int(np.sum(lower*upper%prime,dtype=np.int64))%prime


def scalar_value(colour, tensors_data, prime):
    _,_,lower,upper=tensors_data
    matrices=on_traceless(colour,tensors_data,prime)
    value=upper
    for _ in range(3):
        value=alternate(matrices,value,prime)
    return pairing(lower,value,prime)


def compute(work, prime, points, seed):
    tensor=tensors(prime)
    letters=[(0,0,(i,)) for i in range(3)]
    mu=(0,0,3,3,3)
    images=[m.q(e) for _,block in sorted(m.words_by_shape(8,mu).items()) for e in m.singlets(block)]
    words,rows=m.flatten(images)
    assert all(l in letters for w in words for trace in w for l in trace)
    lib=m.library(work,27,prime)
    f=frame(letters,points,2,0,prime,seed)
    raw=m.values(words,f,lib)
    coeff=np.zeros((len(rows),len(words)),dtype=np.int64)
    for j,row in enumerate(rows):
        for k,v in row:
            coeff[j,k]=v%prime
    assert coeff.shape[1]*(prime-1)**2 < 2**63
    images_at_points=coeff@raw%prime
    polys=[]
    for i in range(points):
        colour=[f.L[l][1].reshape(-1,27,27)[i] for l in letters]
        polys.append(scalar_value(colour,tensor,prime))
    extended=np.vstack((images_at_points,polys))
    minor=m.minor(extended,prime)
    boundary=m.minor(images_at_points,prime)
    return dict(prime=prime,points=points,seed=seed,source_dimension=16,target_dimension=11,
                boundary_minor=boundary,augmented_minor=minor,polynomial_evaluations=polys,
                augmented_evaluations=extended.tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--prime",type=int,default=1000003)
    p.add_argument("--points",type=int,default=16)
    p.add_argument("--seed",type=int,default=614709)
    a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    result=compute(a.work,a.prime,a.points,a.seed)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print("Boundary and augmented ranks",result["boundary_minor"]["rank"],result["augmented_minor"]["rank"])
