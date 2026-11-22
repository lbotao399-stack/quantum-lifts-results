"""Test the complete one-pair terminal highest weight on seven and eight Cartan planes."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations,permutations,product
from math import factorial,prod
from pathlib import Path
from random import Random
import json,sys

from mixed_bmn_d1 import plane_log_kernel,make_class,matmul
from mixed_ward import normalize

def add(x,y):return tuple(a+b for a,b in zip(x,y))
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def sgn(w):return (-1)**sum(x>y for i,x in enumerate(w) for y in w[i+1:])

@lru_cache(None)
def unit(i):return tuple(int(j==i) for j in range(4))

def multiply(a,b):
    row=Counter()
    for w,c in a.items():
        for v,d in b.items():row[add(w,v)]+=c*d
    return row

def dotq(x,y):
    return sum(xi*yi for xi,yi in zip(x,(y[0],y[2],y[1],y[3])))

def quadratic(A):
    q=Counter()
    for i,j in product(range(4),repeat=2):q[add(unit(i),unit(j))]+=dotq(A[i],A[j])
    return q  # twice the quadratic on each plane

def build_probes(rank,fl):
    from flint import nmod_mat
    p=1000003;rng=Random(268511)
    assert fl[0]>=fl[1]>=fl[2] and fl[1]==fl[2]
    eps=fl[2];extra=fl[0]-eps;n=sum(fl)+2
    assert n%2==0 and 2*(rank+2)==n
    triples=list(permutations(range(3)));perms=list(product(triples,repeat=eps))
    signs=[prod(sgn(t) for t in v) for v in perms]
    rows=[];columns={};rank_now=0
    from orthogonal_glueball_endpoints import terminal_weight
    dim=0
    for spin in (0,1):
        for perm in permutations(range(3)):
            w=(1+spin,1-spin)+tuple(fl[i]+perm[i]-i for i in range(3))
            dim+=(-1)**spin*sgn(perm)*sum(c for k,c in terminal_weight(w).items() if k>rank)
    order=[0]*(n//2-2)+[1]
    for attempt in range(600):
        if rank_now==dim:break
        shape=(2,)*(n//2) if attempt%2==0 else (2,)*(n//2-2)+(4,)
        places=list(range(n));rng.shuffle(places)
        site=[i for i,k in enumerate(shape) for _ in range(k)]
        grouped=[[site[j] for j in places[3*i:3*i+3]] for i in range(eps)]
        if any(len(set(t))<3 for t in grouped):continue
        initials=[[0]*4 for k in shape]
        for j in places[3*eps:3*eps+extra]:initials[site[j]][0]+=1
        for j in places[-2:]:initials[site[j]][3]+=1
        row=Counter()
        for perm,sign in zip(perms,signs):
            counts=[u.copy() for u in initials]
            for block,where in zip(perm,grouped):
                for flav,i in zip(block,where):counts[i][flav]+=1
            sites=tuple(sorted(map(tuple,counts),key=lambda w:(sum(w),w)))
            row[sites]+=sign*prod(factorial(k) for w in sites for k in w)
        row={s:c for s,c in row.items() if c}
        if not row:continue
        for w in row:
            if w not in columns:columns[w]=len(columns)
        mat=nmod_mat([[q.get(w,0) for w in columns] for q in rows+[row]],p)
        if mat.rank()>rank_now:rows.append(row);rank_now+=1
    assert rank_now==dim,(rank_now,dim)
    return rows

def finite_functions(seed):
    rng=Random(225677+151757*seed)
    def axis():return [rng.randint(-3,3) for _ in range(4)]
    def sl2():
        a,b,c=axis()[:3];return (a,b,c,-a)
    def form():return {unit(i):n for i,n in enumerate(axis())}
    A=[sl2() for i in range(4)];q=quadratic(A);B=multiply(form(),form());C=multiply(form(),form())
    D=multiply(multiply(form(),form()),form());M=sl2();P=multiply(D,{unit(i):dotq(A[i],M) for i in range(4)})
    BC=multiply(B,C);Bq=multiply(B,q)
    def row(w):
        # q was doubled above. The relative coefficient for two marked paths
        # on one plane and on two planes includes this factor of two.
        if sum(w[-1])==4:
            common=prod(q.get(t,0) for t in w[:-1])
            return (0,0,4*BC.get(w[-1],0)*common,P.get(w[-1],0)*common,Bq.get(w[-1],0)*common)
        a=prod(q.get(t,0) for t in w)
        b=sum(B.get(w[i],0)*prod(q.get(t,0) for j,t in enumerate(w) if i!=j) for i in range(len(w)))
        bc=sum(B.get(u,0)*C.get(v,0)*prod(q.get(t,0) for k,t in enumerate(w) if k not in (i,j))
               for i,u in enumerate(w) for j,v in enumerate(w) if i!=j)
        # Scaling each entire recipe by a constant gives common entries as below:
        # two separate inserts: 4/2^(k-2) versus 2/2^(k-2) on one site.
        return (a,b,8*bc,0,b)
    return row

def bosonic_score(w,score):
    z=[i for i,t in enumerate(w) for _ in range(t[3])]
    assert len(z)==2
    result=0
    for ip,im in set(permutations(z)):
        physical=tuple((int(i==ip),int(i==im))+t[:3] for i,t in enumerate(w))
        s,key=normalize(physical)
        result+=s*(-1)**(im<ip)*score(key)
    return result

def run(a):
    sys.path.insert(0,str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from flint import fmpq_mat
    probes=build_probes(a.rank,tuple(a.flavour))
    fields=sorted({w for r in probes for w in r})
    raw=[];first=[]
    for i in range(a.count):
        op=finite_functions(i);values={w:op(w) for w in fields}
        for j in range(5):
            raw.append([sum(c*values[w][j] for w,c in r.items()) for r in probes])
    print('probes, fields, finite sampled rank',len(probes),len(fields),fmpq_mat(raw).rank(),flush=True)
    log=plane_log_kernel(2*a.rank+1)
    for line in a.parameters.read_text().splitlines()[:4]:
        i,*v=map(int,line.split());A=[(h,e,f,-h) for h,e,f in (v[k:k+3] for k in (0,3,6))]
        score=make_class(A,v[9:12],v[12:15],log)
        values={w:bosonic_score(w,score) for w in fields}
        first.append([str(sum(c*values[w] for w,c in r.items())) for r in probes])
    mat=fmpq_mat(raw+first)
    print('rank with first actions',mat.rank(),flush=True)
    a.out.write_text(json.dumps(dict(physical_flavour=a.flavour,rank=a.rank,probes=len(probes),fields=len(fields),
                  relation_rank=fmpq_mat(raw).rank(),augmented_rank=mat.rank(),finite=raw,first=first),separators=(',',':'))+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rank',type=int,default=6)
    p.add_argument('--flavour',nargs=3,type=int,default=[6,4,4]);p.add_argument('--count',type=int,default=15)
    for k in ('repo','out','parameters'):p.add_argument('--'+k,type=Path,required=True)
    run(p.parse_args())
