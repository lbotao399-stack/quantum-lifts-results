"""The classical map onto purely scalar E6 highest weights in its faithful 27."""
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import sys

import numpy as np
from flint import nmod_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from bmn_primary_index import cochains
import mixed_incoming as m
from f4_colour import basis, jordan_colours, jordan_multiplication
from unitary_trace_weights import highest_weights


@lru_cache(None)
def generators():
    mats = np.concatenate([basis()[0],jordan_colours()])
    assert nmod_mat(mats.reshape(78,-1).tolist(),1000003).rank()==78
    metric = np.array([1]*3+[2]*24,dtype=np.int64)
    unit = np.array([1]*3+[0]*24,dtype=np.int64)
    g = np.diag(metric)
    cubic = np.einsum("i,j,k->ijk",unit,unit,unit)
    linear = np.einsum("i,jk->ijk",unit,g)
    cubic -= linear+linear.transpose(1,2,0)+linear.transpose(2,0,1)
    cubic += jordan_multiplication().transpose(0,2,1)*metric
    assert np.array_equal(cubic,cubic.transpose(1,2,0)) and np.array_equal(cubic,cubic.swapaxes(0,1))
    for mat in mats:
        deriv = np.einsum("ai,ajk->ijk",mat,cubic)
        assert not (deriv+deriv.transpose(1,2,0)+deriv.transpose(2,0,1)).any()
    return mats


def frame(letters, points, prime, seed):
    mats=generators()
    rng=np.random.default_rng(seed)
    evals=m.Frame(27,points,1,[0],letters,seed=1,p=prime)
    for l in sorted(letters):
        parity,template=evals.L[l]
        coordinates=rng.integers(0,prime,size=template.shape[:-2]+(78,),dtype=np.int64)
        evals.L[l]=parity,(coordinates@mats.reshape(78,-1)).reshape(template.shape)%prime
    return evals


def trace_boundaries(weight, parity, evals, lib, prime):
    parents = highest_weights((0,0,*weight), sum(weight)-1, parity)
    images = []
    for source in parents.get(1, ()):
        image = Counter()
        for word, coefficient in source.items():
            for target, n in m.fw.Q_multi(word).items():
                image[target] += coefficient*n
        images.append({w:c for w,c in image.items() if c})
    words, rows = m.flatten(images)
    matrix = np.zeros((len(rows), evals.ncols), dtype=np.int64)
    if words:
        values = m.values(words, evals, lib)
        for j, row in enumerate(rows):
            for k, coefficient in row:
                matrix[j] = (matrix[j]+(coefficient%prime)*values[k])%prime
    return matrix


def compute(weights, points, prime, seed, work):
    work.mkdir(parents=True,exist_ok=True)
    lib=m.library(work,27,prime)
    letters=[(0,0,(i,)) for i in range(3)]
    evals=frame(letters,points,prime,seed)
    dimensions=cochains(max(map(sum,weights)),"E6")[1]
    records=[]
    for w in weights:
        d=sum(w)
        for i,parity in enumerate((1,-1)):
            dim=dimensions[d].get(w,[0,0])[i]
            if dim==0:
                continue
            matrix=trace_boundaries(w,parity,evals,lib,prime)
            witness=m.minor(matrix,prime)
            assert witness["rank"]<=dim
            records.append(dict(weight=w,length=d,outer=parity,dimension=dim,parents=len(matrix),
                                cohomology_upper=dim-witness["rank"],boundary_minor=witness,
                                minor_matrix=matrix[np.ix_(witness["rows"],witness["columns"])].tolist()))
            print(w,parity,dim,len(matrix),witness["rank"], flush=True)
    return dict(prime=prime,seed=seed,points=points,verified_cubic_generators=78,results=records)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--prime",type=int,choices=(1000003,1000033),default=1000003)
    p.add_argument("--seed",type=int,default=334057)
    p.add_argument("--points",type=int,default=22)
    p.add_argument("--weight",action="append",nargs=3,type=int,required=True)
    a=p.parse_args()
    data=compute([tuple(w) for w in a.weight],a.points,a.prime,a.seed,a.work)
    a.out.write_text(json.dumps(data,separators=(",",":"))+"\n")
