"""Verify the fixed 703-minors, the character dimensions and the Cartan complement."""
from argparse import ArgumentParser
from collections import Counter
from functools import cache
from hashlib import sha256
from pathlib import Path
import gzip
import json

import numpy as np

HERE=Path(__file__).resolve().parent


def det_mod(array,p):
    m=np.array(array,dtype=np.int64,copy=True)
    det=1
    for i in range(len(m)):
        indices=np.flatnonzero(m[i:,i])
        if not len(indices):return 0
        j=i+int(indices[0])
        if j!=i:
            m[[i,j]]=m[[j,i]]
            det=-det
        value=int(m[i,i]);det=det*value%p
        co=m[i+1:,i]*pow(value,-1,p)%p
        m[i+1:,i:]=(m[i+1:,i:]-co[:,None]*m[i,i:])%p
    return det


def integer_det(m):
    from itertools import permutations
    from math import prod
    return sum((-1)**sum(a>b for i,a in enumerate(perm) for b in perm[i+1:])
               *prod(int(m[i,j]) for i,j in enumerate(perm)) for perm in permutations(range(len(m))))


def characters_and_cartan():
    from b2_schur_complex import multiplicity
    from f4_singlet_characters import decomposition
    @cache
    def singlet(cycles):
        a,b=[],[]
        for t in sorted(cycles,reverse=True):
            (a if sum(a)<=sum(b) else b).append(t)
        left,right=decomposition("F4",tuple(a)),decomposition("F4",tuple(b))
        return sum(v*right.get(k,0) for k,v in left.items())
    dims=[multiplicity((4,4,3,3),j,singlet) for j in (12,13,14)]
    assert dims==[10007,3720,704]
    raw=json.loads(gzip.decompress((HERE/"representation26.json.gz").read_bytes()))["matrices"]
    diag=np.zeros((4,26),dtype=np.int64)
    for i,row in enumerate(raw[:4]):
        for j,k,v in row:
            assert j==k
            diag[i,j]=v
    gram=diag@diag.T
    dets=[integer_det(gram[:n,:n]) for n in (3,4)]
    assert dets[0]*dets[1]!=0
    return dims,dets


def read_matrix(p):
    info=json.loads((HERE/f"selected-{p}.json").read_bytes())
    raw=gzip.decompress((HERE/f"selected-{p}.u16.gz").read_bytes())
    assert sha256(raw).hexdigest()==info["residue_u16le_sha256"]
    assert len(raw)==2*703**2 and len(set(info["points"]))==703
    assert len({tuple(v) for v in info["rows"]})==703
    matrix=np.frombuffer(raw,dtype="<u2").reshape(703,703)
    assert int(matrix.max())<p
    return info,matrix


def verify(a):
    dimensions,gram=characters_and_cartan()
    records={}
    baseline=None
    for prime in (1009,1013):
        info,matrix=read_matrix(prime)
        identity=info["rows"],info["points"]
        if baseline is None:baseline=identity
        assert baseline==identity
        det=det_mod(matrix,prime)
        assert det==info["determinant"]=={1009:1005,1013:734}[prime]
        counts=dict(Counter(stage for stage,i in info["rows"]))
        assert counts=={"reference":542,"four":158,"nested":3}
        records[str(prime)]=det
        if a.samples or a.full:
            from rebuild_classical import regenerate
            positions=range(703) if a.full else a.samples
            for pos,values in regenerate(prime,info,positions):
                assert np.array_equal(matrix[:,pos],values), (prime,pos)
    return dict(dimensions=dimensions,cartan_gram_determinants=gram,
                cartan_value=gram[0]*gram[1],rows_by_family=counts,minor_determinants=records,
                classical_rank=703,classical_kernel_dimension=dimensions[1]-703,
                terminal_cokernel_dimension=dimensions[2]-703)


if __name__=="__main__":
    p=ArgumentParser(description=__doc__)
    mode=p.add_mutually_exclusive_group()
    mode.add_argument("--samples",type=int,nargs="+");mode.add_argument("--full",action="store_true")
    a=p.parse_args();assert __debug__
    print(json.dumps(verify(a),indent=2))
