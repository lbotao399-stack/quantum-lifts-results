"""Two exceptional invariant contractions of open words in the minimal F4 representation."""
from fractions import Fraction
from pathlib import Path
import gzip,json
import numpy as np


def projectors(folder, rep, prime):
    # End(26) has five nonisomorphic F4 summands with these quadratic Casimirs.
    # All denominators in the native inverse metric divide six or a small integer.
    data=json.loads(gzip.decompress((folder/'inputs.json.gz').read_bytes()))
    ts=rep.astype(np.int64)
    dim=26
    ident=np.eye(dim,dtype=np.int64)
    basic=np.zeros_like(ident)
    L=np.zeros((dim**2,dim**2),dtype=np.int64)
    for i,j,entry in data['inverse']:
        q=Fraction(entry)
        v=q.numerator*pow(q.denominator,-1,prime)%prime
        basic=(basic+v*(ts[i]@ts[j]))%prime
        # Row-major flattening of T_i X T_j.
        L=(L-2*v*np.kron(ts[i],ts[j].T))%prime
    c=int(basic[0,0]);assert np.array_equal(basic,c*ident)
    assert c==12
    L=(L+2*c*np.eye(dim**2,dtype=np.int64))%prime
    eig=(0,12,18,24,26)
    expected=(1,26,52,273,324)
    def make(k):
        out=np.eye(dim**2,dtype=np.int64)
        for v in eig:
            if v!=k:
                out=(L@out-v*out)*pow(k-v,-1,prime)%prime
        return out
    result={}
    for k,requested,transpose_sign in ((12,26,1),(18,52,-1)):
        proj=make(k)
        assert (L@proj-k*proj).max()>=0
        assert not np.any((L@proj-k*proj)%prime)
        selected=[];basis={}
        for j,col in enumerate(proj.T):
            x=col.copy()
            for i in sorted(basis):
                if x[i]:x=(x-x[i]*basis[i])%prime
            nz=np.flatnonzero(x)
            if len(nz):
                i=int(nz[0]);basis[i]=x*pow(int(x[i]),-1,prime)%prime;selected.append(j)
        assert len(basis)==requested
        J=proj[:,selected].T.reshape(-1,dim,dim)
        # Transposition is in the invariant bilinear form on 26.
        minimal=json.loads(gzip.decompress((folder/'representation26.json.gz').read_bytes()))
        h=np.zeros((dim,dim),dtype=np.int64)
        for i,j,t in minimal['form']:
            q=Fraction(t)
            h[i,j]=q.numerator*pow(q.denominator,-1,prime)%prime
        assert not np.any((h@J-transpose_sign*J.transpose(0,2,1)@h)%prime)
        vectors=J.transpose(0,2,1).reshape(requested,-1)
        gram=vectors@J.reshape(requested,-1).T%prime
        A=np.column_stack((gram,np.eye(requested,dtype=np.int64)))
        for i in range(requested):
            j=i+int(np.flatnonzero(A[i:,i])[0]);A[[i,j]]=A[[j,i]]
            A[i]=A[i]*pow(int(A[i,i]),-1,prime)%prime
            co=A[:,i].copy();co[i]=0;A=(A-co[:,None]*A[i])%prime
        inv=A[:,requested:]
        assert np.array_equal(gram@inv%prime,np.eye(requested,dtype=np.int64))
        result[k]=vectors,inv
    return result
