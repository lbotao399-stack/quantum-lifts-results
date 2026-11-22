"""Evaluate the six odd commutant substitutions and their return on first-jet traces."""
import argparse
from collections import Counter
from functools import lru_cache
from math import factorial,prod
import json
from pathlib import Path
import sys


def normalize(word):
    def key(x): return sum(x),x
    odds=[key(w) for w in word if (w[0]+w[1])%2]
    if len(set(odds))<len(odds):return 0,()
    s=(-1)**sum(a>b for i,a in enumerate(odds) for b in odds[i+1:])
    return s,tuple(sorted(word,key=key))


@lru_cache(None)
def operation(w,spin,flavour,reverse):
    result=Counter()
    prefix=0
    for i,t in enumerate(w):
        a=list(t)
        c=(a[flavour+2] if not a[spin] else 0) if reverse else a[spin]
        if c:
            a[spin]+=1 if reverse else -1
            a[flavour+2]+=-1 if reverse else 1
            s,v=normalize(w[:i]+(tuple(a),)+w[i+1:])
            if s:
                result[v]+=s*c*(-1)**(prefix+sum(t[:spin]))
        prefix+=t[0]+t[1]
    return tuple((v,c) for v,c in result.items() if c)


def transform(expr,spin,flavour,reverse):
    ret=Counter()
    for w,c in expr.items():
        for v,x in operation(w,spin,flavour,reverse):
            ret[v]+=c*x
    return {w:c for w,c in ret.items() if c}


def return_map(basis,steps=None):
    from flint import fmpq_mat
    if steps is None:
        steps=tuple((s,f) for s in range(2) for f in range(3))
    matrix=[]
    for w in basis:
        expr={w:1}
        for s,f in steps: expr=transform(expr,s,f,False)
        for s,f in reversed(steps): expr=transform(expr,s,f,True)
        assert set(expr)<=set(basis)
        matrix.append([expr.get(v,0) for v in basis])
    return fmpq_mat(matrix)


def trace_to_planes(basis,n):
    from flint import fmpq_mat
    from orthogonal_glueball_cartan_basis import groupings
    plane_list=[w for w in basis if len(w)<=n]
    index={w:i for i,w in enumerate(plane_list)}
    rows=[]
    for w in basis:
        row=[0]*len(plane_list)
        factor=2**len(w)*(-1)**(sum(map(sum,w))//2)
        for groups in groupings(len(w)):
            if len(groups)>n:continue
            blocks=[tuple(sum(w[i][j] for i in g) for j in range(5)) for g in groups]
            if any(max(b[:2])>1 for b in blocks):continue
            s,term=normalize(tuple(blocks))
            if not s:continue
            order=sorted(range(len(blocks)),key=lambda i:(sum(blocks[i]),blocks[i]))
            labels={k:j for j,g in enumerate(order) for k in groups[g]}
            seq=[2*labels[i]+k for i,t in enumerate(w) for k in range(2) if t[k]]
            sign=(-1)**sum(t>u for i,t in enumerate(seq) for u in seq[i+1:])
            row[index[term]]+=factor*sign*prod(factorial(k) for k in Counter(blocks).values())
        rows.append(row)
    return fmpq_mat(rows)


def kernel_rows(matrix):
    from flint import fmpq_mat
    red,rank=matrix.transpose().rref()
    pivots=[next(j for j in range(red.ncols()) if red[i,j]) for i in range(rank)]
    rows=[]
    for k in range(red.ncols()):
        if k in pivots: continue
        v=[0]*red.ncols()
        v[k]=1
        for i,j in enumerate(pivots):v[j]=-red[i,k]
        rows.append(v)
    mat=fmpq_mat(rows) if rows else fmpq_mat(0,matrix.nrows())
    assert not any((mat*matrix).entries())
    return mat


def format_matrix(m):return [[str(v) for v in row] for row in m.tolist()]


def run(case):
    from flint import fmpq_mat
    n=case['colours']//2
    basis=[tuple(map(tuple,w)) for w in case['target']]
    res=trace_to_planes(basis,n)
    kernel=kernel_rows(res)
    ward=return_map(basis)
    # In the pure powers the independent native six-step calculation gives these exact values.
    pure=[((1,1,0,0,0),)*p for p in (5,6)]
    assert [int(return_map([w])[0,0]) for w in pure]==[7200,25200]
    boundaries=fmpq_mat(case['image_basis']) if case.get('image_basis') else fmpq_mat(0,len(basis))
    assert not any((boundaries*res).entries()) and not any((kernel*ward*res).entries())
    assert setrow(boundaries,boundaries*ward).rank()==boundaries.rank()
    partial={}
    for spin in (0,1):
        for mask in range(8):
            steps=tuple((spin,f) for f in range(3))+tuple((1-spin,f) for f in range(3) if mask&(1<<f))
            ret=return_map(basis,steps)
            assert not any((kernel*ret*res).entries())
            assert setrow(boundaries,boundaries*ret).rank()==boundaries.rank()
            partial[f'{spin}:{mask}']=setrow(boundaries,kernel*ret).rank()-boundaries.rank()
    return dict(colours=case['colours'],weight=case['weight'],formal=len(basis),cartan_rank=res.rank(),
                cartan= format_matrix(res),cartan_kernel=format_matrix(kernel),ward=format_matrix(ward),
                d1_rank=boundaries.rank(),return_image_on_cartan_kernel_mod_d1=
                setrow(boundaries,kernel*ward).rank()-boundaries.rank(),
                one_spin_three_returns=partial,
                remaining_cartan_kernel=kernel.rank()-boundaries.rank())


def setrow(a,b):
    from flint import fmpq_mat
    rows=a.tolist()+b.tolist()
    return fmpq_mat(rows) if rows else fmpq_mat(0,a.ncols())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',required=True,type=Path)
    p.add_argument('--case',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    assert __debug__
    sys.path.insert(0,str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    d=run(json.loads(a.case.read_text()))
    a.out.write_text(json.dumps(d,indent=2)+'\n')
    print({k:v for k,v in d.items() if k not in ('cartan','cartan_kernel','ward')})
