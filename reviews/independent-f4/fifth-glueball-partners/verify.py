"""Regenerate the two rational rank bounds using integer arithmetic."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from functools import cache
import gzip
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix

import exceptional_colour
import recipes
from b2_schur_complex import multiplicity
from f4_singlet_characters import decomposition

HERE = Path(__file__).resolve().parent


@cache
def singlet(cycles):
    a,b=[],[]
    for n in sorted(cycles,reverse=True):
        (a if sum(a)<=sum(b) else b).append(n)
    x=decomposition("F4",tuple(a))
    y=decomposition("F4",tuple(b))
    return sum(v*y.get(k,0) for k,v in x.items())


def ordinary_values(frame, coeff, batch, p):
    matrix, groups, factors, atoms = batch
    values=np.zeros(atoms,dtype=np.int64)
    for i,rep in enumerate(frame):
        size=rep.shape[1]
        letters=(coeff@rep.reshape(52,-1)).reshape(10,size,size)%p
        local={k:v for k,v in groups.items() if k[0]==i}
        prefixes=defaultdict(set)
        for (_,_,length),(_,codes) in local.items():
            m=length//2
            prefixes[m].update(map(int,codes//4**(length-m)))
            prefixes[length-m].update(map(int,codes%4**(length-m)))
        for n in range(max(prefixes,default=0),0,-1):
            prefixes[n-1].update(k//4 for k in prefixes[n])
        tables={0:np.eye(size,dtype=np.int64)[None]}
        indices={0:{0:0}}
        for n in range(1,max(prefixes,default=0)+1):
            labels=sorted(prefixes[n])
            ix=np.asarray([indices[n-1][k//4] for k in labels])
            matrices=tables[n-1][ix] @ letters[np.asarray(labels)%4]
            tables[n]=matrices%p
            indices[n]={k:j for j,k in enumerate(labels)}
        for (_,pair,length),(positions,codes) in local.items():
            m=length//2
            lo=np.asarray([indices[m][int(k)//4**(length-m)] for k in codes])
            hi=np.asarray([indices[length-m][int(k)%4**(length-m)] for k in codes])
            left=tables[m] if pair==0 else (letters[3+pair]@tables[m])%p
            right=tables[length-m]
            for j in range(0,len(lo),80):
                s=slice(j,j+80)
                vals=np.einsum("nij,nji->n",left[lo[s]],right[hi[s]],dtype=np.int64)
                values[positions[s]]=vals%p
    result=np.zeros(matrix.shape[1],dtype=np.int64)
    for ids,words in factors.values():
        v=np.ones(len(ids),dtype=np.int64)
        for factor in words.T:
            v=v*values[factor]%p
        result[ids]=v
    return np.asarray(matrix@result)%p


def projected_values(rep,coeff,batch,dual,inv,p):
    matrix,types,layers,groups,opened,closed=batch
    assert set(types)=={12} and not groups and closed==1
    letters=(coeff@rep.reshape(52,-1)).reshape(10,26,26)%p
    projected=np.zeros((opened,26),dtype=np.int64)
    previous=np.eye(26,dtype=np.int64)[None]
    for prefix,last,where,dest in layers:
        current=(previous[prefix]@letters[last])%p
        if len(where):
            projected[dest]=current[where].reshape(-1,676)@dual.T%p
        previous=current
    ids,a,b,c=types[12]
    assert not np.any(c)
    v=np.zeros(matrix.shape[1],dtype=np.int64)
    v[ids]=np.einsum("ij,ij->i",(projected@inv%p)[a],projected[b])%p
    return np.asarray(matrix@v)%p


def signed(order):
    return (-1)**sum(a>b for n,a in enumerate(order) for b in order[n+1:])


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


def nested_rows(rep,coeff,ends,p):
    low,high,transfer=ends
    x=(coeff[:4]@rep.reshape(52,-1)).reshape(4,26,26)%p
    brackets={(i,j):(x[i]@x[j]-x[j]@x[i])%p for i,j in combinations(range(4),2)}
    def end(v,is_lower):
        xx=[csr_matrix(m.T if is_lower else m) for m in x]
        out=np.zeros_like(v)
        for (i,j),m in brackets.items():
            k,l=(q for q in range(4) if q not in (i,j))
            t=act(csr_matrix(m.T if is_lower else m),v,0,p)
            a=act(xx[l],act(xx[k],t,1,p),2,p)
            b=act(xx[k],act(xx[l],t,1,p),2,p)
            out+=signed((i,j,k,l))*(a-b)
        return out%p
    left=end(low,True).ravel()
    right=end(high,False)
    xx=[csr_matrix(m) for m in x]
    ordinary=0
    exceptional=0
    for i,j,k in permutations(range(3)):
        t=act(xx[i],right,0,p)
        u=act(transfer,t,(0,1),p)
        for which,v in enumerate((t,u)):
            v=act(xx[k],act(xx[j],v,1,p),2,p)
            value=signed((i,j,k))*int(left@v.ravel())
            if which==0:
                ordinary+=value
            else:
                exceptional+=value
    return [ordinary%p,exceptional%p]


def bareiss_residue(matrix,p):
    a=[[int(v) for v in row] for row in matrix]
    last=1
    sign=1
    for i in range(len(a)-1):
        if a[i][i]==0:
            j=next(k for k in range(i+1,len(a)) if a[k][i])
            a[i],a[j]=a[j],a[i]
            sign=-sign
        pivot=a[i][i]
        for j in range(i+1,len(a)):
            lead=a[j][i]
            for k in range(i+1,len(a)):
                value=a[j][k]*pivot-lead*a[i][k]
                assert value%last==0
                a[j][k]=value//last
        last=pivot
    return sign*a[-1][-1]%p


def minor(mode,p,frame,dual,inv,ends):
    seeds=json.loads((HERE/"recipes.json").read_bytes())[mode]
    batch=recipes.inputs(mode,seeds["trace_recipes"],True)
    extra=recipes.prepare(mode,seeds["exceptional_recipes"])
    n=102 if mode=="in" else 36
    cols=[]
    for i in range(n):
        rng=np.random.Generator(np.random.PCG64(763931+307*i))
        coeff=rng.integers(-4,5,size=(10,52),dtype=np.int64)
        last=nested_rows(frame[0],coeff,ends,p) if mode=="out" else ()
        cols.append(np.concatenate((ordinary_values(frame,coeff,batch,p),
                                    projected_values(frame[0],coeff,extra,dual,inv,p),
                                    np.asarray(last,dtype=np.int64))))
        if (i+1)%12==0:
            print(mode,p,"points",i+1,file=sys.stderr,flush=True)
    a=np.stack(cols,axis=1)
    assert a.shape==(n,n)
    det=bareiss_residue(a,p)
    rs,cs,pivot_product=recipes.pivots(a,p)
    assert len(rs)==n and (signed(cs)*pivot_product)%p==det and det
    from hashlib import sha256
    digest=sha256(a.astype("<u2").tobytes()).hexdigest()
    if p==1009:
        expected={"in":(313,"3b7d4018518c7028c50f9d129fd7c03ece24ff0aa2e6df7fe6aa85f664c3d7a6"),
                  "out":(987,"55d6bbc097c8d8b8915bf37ff517b76bff0408f84a3fc6e95bafaf3bc3796342")}
        assert (det,digest)==expected[mode],(mode,det,digest)
    return dict(order=n,determinant=det,entries_sha256=digest)


def run(p,mode):
    assert p in (1009,1013)
    # The widest dot product contains 26^4 entries smaller than p in magnitude.
    assert 26**4*(p-1)**2<2**63
    dims=[multiplicity((4,4,3),y,singlet) for y in (9,10,11)]
    assert dims==[317,138,36]
    frame=recipes.matrices(HERE)
    dual,inv=exceptional_colour.projectors(HERE,frame[0],p)[12]
    ends=invariant_ends(frame[0],dual,inv,p)
    selected=("in","out") if mode=="both" else (mode,)
    return dict(prime=p,dimensions=dims,
                minors={k:minor(k,p,frame,dual,inv,ends) for k in selected})


if __name__=="__main__":
    a=ArgumentParser(description=__doc__)
    a.add_argument("--prime",type=int,choices=(1009,1013),required=True)
    a.add_argument("--mode",choices=("both","in","out"),default="both")
    if not __debug__:
        a.error("run with assertions enabled")
    opts=a.parse_args()
    print(json.dumps(run(opts.prime,opts.mode),indent=2))
