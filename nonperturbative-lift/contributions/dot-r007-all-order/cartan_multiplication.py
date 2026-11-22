"""Exact small Cartan map for the two R16_[0,1] coefficients in dMK."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

# A component letter: dotted exponents a,b, theta mask, independent Cartan color 0 or 1.
def odd(x):return (x[2].bit_count()+1)%2
def canon(w):
    o=[x for x in w if odd(x)]
    if len(set(o))<len(o): return (),0
    return tuple(sorted(w)),(-1)**sum(x>y for i,x in enumerate(o) for y in o[i+1:])
def add(*terms):
    R=defaultdict(F)
    for c,A in terms:
        for w,a in A.items():
            if c:R[w]+=a*c
    return {w:c for w,c in R.items() if c}
def mul(A,B):
    R=defaultdict(F)
    for w,v in product(A,B):
        z,s=canon(w+v)
        if s:R[z]+=s*A[w]*B[v]
    return {w:c for w,c in R.items() if c}
def atom(mask,c):
    return {((0,0,mask,c),):F(1)}
def at(mask,row):
    if row<2:return atom(mask,row)
    return add((-1,atom(mask,0)),(-1,atom(mask,1)))
def trace_simple(*masks):
    T=[]
    for a in range(3):
        P={():F(1)}
        for m in masks:P=mul(P,at(m,a))
        T.append((1,P))
    return add(*T)
@lru_cache(None)
def raw(i,k):
    return add(((-1)**k,trace_simple(1<<i,7^(1<<k))))
@lru_cache(None)
def v2(i,j):
    return add((1,raw(i,j)),*[(F(-1,3),raw(k,k)) for k in range(3) if i==j])
@lru_cache(None)
def raw3(i,j,k):
    return add(((-1)**k,trace_simple(1<<i,1<<j,7^(1<<k))))
@lru_cache(None)
def v3(i,j,k):
    return add((1,raw3(i,j,k)),*[(F(-1,4),raw3(j,a,a)) for a in range(3) if i==k],*[(F(-1,4),raw3(i,a,a)) for a in range(3) if j==k])
def eps(a,b,c):
    return (-1)**((a>b)+(a>c)+(b>c)) if len({a,b,c})==3 else 0
coeff=[{} for _ in range(3)]
for a,b,c in product(range(3),repeat=3):
    e=eps(a,b,0)
    if not e:continue
    coeff[b]=add((1,coeff[b]),(F(8,5)*e,mul(v2(c,a),v3(0,0,c))))
    coeff[c]=add((1,coeff[c]),(F(16,15)*e,mul(v2(0,a),v3(0,c,b))))
print('prefactors',[len(Fi) for Fi in coeff],flush=True)
t=(2,2,1,1,0)
Ls=[(a,b,m,c) for a,b,m,c in product(range(3),range(3),range(4),range(2)) if a+b+m]
wt=lambda x:(x[0],x[1],*(int(bool(x[2]&(1<<i))) for i in range(3)))
Lw=list(map(wt,Ls))
@lru_cache(None)
def mns(ia,degree,W):
    if not degree:return ((),) if not any(W) else ()
    if sum(W)<degree:return ()
    ret=[]
    for j in range(ia,len(Ls)):
        rem=tuple(a-b for a,b in zip(W,Lw[j]))
        if min(rem)<0:continue
        for w in mns(j+odd(Ls[j]),degree-1,rem):
            ret.append((Ls[j],)+w)
    return tuple(ret)
def flavor(w,k):
    new=[]
    s=1
    for a,b,m,c in w:
        seq=[(i+k)%3 for i in range(3) if m&(1<<i)]
        s*=(-1)**sum(x>y for j,x in enumerate(seq) for y in seq[j+1:])
        new.append((a,b,sum(1<<i for i in seq),c))
    z,h=canon(new)
    return z,s*h
basis=mns(0,4,t)
print('raw balanced highest-weight monomials before projection',len(basis),flush=True)
def gauge_swap(w,i):
    # Simple Weyl s0 exchanges colors 0,1; s1 exchanges diagonal entries 1,2.
    P={():F(1)}
    for a,b,m,c in w:
        A={((a,b,m,1-c),):F(1)} if i==0 else ({((a,b,m,0),):F(1)} if c==0 else {((a,b,m,0),):F(-1),((a,b,m,1),):F(-1)})
        P=mul(P,A)
    return add((1,P),(-1,{w:F(1)}))
def raising(w,i):
    out=defaultdict(F)
    for j,(a,b,m,c) in enumerate(w):
        if i==2:
            if b==0:continue
            v=(a+1,b-1,m,c);d=a+1
        else:
            if not (m&(1<<(i+1))) or m&(1<<i):continue
            v=(a,b,m^(3<<i),c);d=1
        z,s=canon(w[:j]+(v,)+w[j+1:])
        if s:out[z]+=s*d
    return {w:c for w,c in out.items() if c}
def map_column(w):
    z0,s0=flavor(w,1);z1,s1=flavor(w,2)
    return add((1,mul(coeff[2],{w:F(1)})),(s0,mul(coeff[0],{z0:F(1)})),(s1,mul(coeff[1],{z1:F(1)})))
def mod(v,p):return v.numerator % p * pow(v.denominator,-1,p) % p
def rank(columns,p):
    piv={}
    for P in columns:
        col={w:mod(v,p) for w,v in P.items() if mod(v,p)}
        while col:
            x=min(col);c=col[x]
            if x not in piv:
                iv=pow(c,-1,p);piv[x]={w:v*iv%p for w,v in col.items()}
                break
            for w,v in piv[x].items():
                n=(col.get(w,0)-c*v)%p
                if n:col[w]=n
                else:col.pop(w,None)
    return len(piv)
constraints=[]
M=[]
for w in basis:
    ops=[gauge_swap(w,i) for i in range(2)]+[raising(w,i) for i in range(3)]
    constraints.append({(i,z):c for i,P in enumerate(ops) for z,c in P.items()})
    M.append({(5,z):c for z,c in map_column(w).items()})
for p in (1000003,2147483647):
    r=rank(constraints,p)
    j=rank([A|B for A,B in zip(constraints,M)],p)
    print({'prime':p,'raw':len(basis),'constraint_rank':r,'ambient_hw_gauge_dim':len(basis)-r,'stacked_rank':j,'multiplication_kernel_upper_Q':len(basis)-j},flush=True)
