"""Integer Albert endpoint actions from the fifth-power partner calculation."""
from fractions import Fraction
from pathlib import Path
import gzip
import json
import numpy as np
from scipy.sparse import csr_matrix

HERE=Path(__file__).resolve().parent

def act(a,v,axis,p):
    # Integer CSR multiplication on one or two specified tensor slots.
    axes=(axis,) if isinstance(axis,int) else axis
    front=tuple(range(len(axes)))
    aligned=np.ascontiguousarray(np.moveaxis(v,axes,front))
    size=26**len(axes)
    result=(a@aligned.reshape(size,-1))%p
    return np.moveaxis(result.reshape((26,)*4),front,axes)

def invariant_ends(rep,coords,inv,p):
    h=np.zeros((26,26),dtype=np.int64)
    for i,j,value in json.loads(gzip.decompress((HERE/"representation26.json.gz").read_bytes()))["form"]:
        q=Fraction(value)
        h[i,j]=q.numerator*pow(q.denominator,-1,p)%p
    hi=h.copy()
    hi[:2,:2]=np.array([[2,1],[1,2]])*pow(3,-1,p)%p
    assert np.array_equal(h@hi%p,np.eye(26,dtype=np.int64))
    jordan=coords.reshape(-1,26,26).transpose(0,2,1)
    lo=(h@jordan%p).reshape(26,-1)
    up=(jordan@hi%p).reshape(26,-1)
    low=((lo.T@inv%p@lo)%p).reshape((26,)*4)
    high=((up.T@inv%p@up)%p).reshape((26,)*4).transpose(0,2,1,3)
    raw=jordan.reshape(26,-1)
    transfer=((raw.T@inv%p@raw)%p).reshape((26,)*4).transpose(0,2,1,3)
    transfer=csr_matrix(transfer.reshape(676,676))
    eye=np.eye(26,dtype=np.int64)
    for i in (27,51,10,24,28,4,43,30):
        for v,a in ((low,rep[i].T),(high,rep[i])):
            assert not np.any(sum(act(csr_matrix(a),v,k,p) for k in range(4))%p)
        adj=csr_matrix(np.kron(rep[i],eye)+np.kron(eye,rep[i]))
        assert not np.any((adj@transfer-transfer@adj).data%p)
    return low,high,transfer
