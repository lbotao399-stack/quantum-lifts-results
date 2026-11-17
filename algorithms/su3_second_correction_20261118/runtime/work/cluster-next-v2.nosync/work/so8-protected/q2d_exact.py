"""Exact two-loop trace action with input jets, in the pinned upstream convention.

Theta/color rules extend davidsd/fortuity-ai commit
4256f207986029540060fb44dd494746828f65ce, q2bmn.py.
Input-jet moments use n4_supercharges.tex, eqs. kernel/currents2/lowmeasures:
integral dt1 dt2 dt3 dt4 t1 t2 t4, t3+t4<=1.
Its zero-jet moment is 1/24. Relative moments retain upstream K2=-1/6.
SO traces use both terms of the same completeness relation as q1.py.
"""
import sys,pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent.parent/'upstream'))
from q2bmn import triple_rule,_reorder_sign
from functools import lru_cache
from itertools import product,combinations
from fractions import Fraction
from math import factorial
from collections import defaultdict
from fortuity.full.letters import is_ferm
from fortuity.full.words import canon_word,canon_multi

ONE={(0,0,0,0):1}


def add(*terms):
    out=defaultdict(int)
    for scalar,p in terms:
        for m,c in p.items():out[m]+=scalar*c
    return {m:c for m,c in out.items() if c}


def mul(a,b):
    out=defaultdict(int)
    for m,c in a.items():
        for n,d in b.items():out[tuple(x+y for x,y in zip(m,n))]+=c*d
    return {m:c for m,c in out.items() if c}


T=[{tuple(int(i==j) for i in range(4)):1} for j in range(4)]
NT=[add((1,ONE),(-1,t)) for t in T]
A=(add((1,NT[0]),(1,mul(mul(T[0],NT[1]),NT[2]))),
   mul(mul(T[0],T[1]),NT[2]),mul(T[0],T[2]))
B=(mul(NT[1],T[3]),mul(T[1],T[3]),NT[3])


@lru_cache(None)
def moment(exponents):
    pol=ONE
    for p,n in zip(A+B,exponents):
        for _ in range(n):pol=mul(pol,p)
    # Independent intervals t1,t2 and triangle t3,t4, including t1*t2*t4.
    value=Fraction(0)
    for (i,j,k,l),c in pol.items():
        value+=Fraction(c*factorial(k)*factorial(l+1),(i+2)*(j+2)*factorial(k+l+3))
    return 24*value


@lru_cache(None)
def jet_terms(jets):
    out=defaultdict(Fraction)
    choices=[list(product(range(a+1),range(b+1))) for a,b in jets]
    for split in product(*choices):
        d1=sum(x[0] for x in split);d2=sum(x[1] for x in split)
        e1=sum(x[0] for x in jets)-d1;e2=sum(x[1] for x in jets)-d2
        den=1
        for (a,b),(r,s) in zip(jets,split):den*=factorial(r)*factorial(s)*factorial(a-r)*factorial(b-s)
        exponents=tuple(r+s for r,s in split)+tuple(a+b-r-s for (a,b),(r,s) in zip(jets,split))
        w=moment(exponents)/den
        for scalar,da,db,ea,eb in ((1,2,0,0,2),(-2,1,1,1,1),(1,0,2,2,0)):
            key=(d1+da,d2+db,e1+ea,e2+eb)
            cf=scalar*w
            for n in key:cf*=factorial(n)
            out[key]+=cf
    return tuple((c,(a,b),(d,e)) for (a,b,d,e),c in out.items() if c)


def contract(necks,tok,group):
    (p,u),(q,v)=[(a,i) for a,w in enumerate(necks) for i,t in enumerate(w) if t==tok]
    if p==q:
        x=necks[p];aa=x[u+1:v];bb=x[v+1:]+x[:u];rest=[w for a,w in enumerate(necks) if a!=p]
        yield 1,rest+[aa,bb]
        if group=='SO_odd':yield -(-1)**len(aa),rest+[aa[::-1]+bb]
    else:
        x,y=necks[p],necks[q];aa=x[u+1:]+x[:u];bb=y[v+1:]+y[:v]
        rest=[w for a,w in enumerate(necks) if a not in (p,q)]
        yield 1,rest+[aa+bb]
        if group=='SO_odd':yield -(-1)**len(bb),rest+[aa+bb[::-1]]


def Q2d_multi(multi,N,group='SU'):
    if group not in ('SU','SO_odd'):raise NotImplementedError(group)
    letters=[L for w in multi for L in w];n=len(letters)
    par=[int(is_ferm(L)) for L in letters];ids=[];k=0
    for w in multi:ids.append(list(range(k,k+len(w))));k+=len(w)
    out=defaultdict(Fraction)
    for tri in combinations(range(n),3):
        rest=[t for t in range(n) if t not in tri]
        for middle in tri:
            left,right=[t for t in tri if t!=middle]
            rule=triple_rule(letters[left][2],letters[middle][2],letters[right][2])
            if not rule:continue
            sf=_reorder_sign(list(range(n)),[left,middle,right]+rest,lambda t:par[t])
            jets=jet_terms(tuple(letters[t][:2] for t in (left,middle,right)))
            ref=[left,right]+rest
            for c0,SD,SE in rule:
                newpar={left:(len(SD)+1)%2,right:(len(SE)+1)%2}
                parity=lambda t:newpar.get(t,par[t])
                for ca,cm,cc in product((0,1),repeat=3):
                    necks=[]
                    for word in ids:
                        w=[]
                        for t in word:
                            if t==left:w+=['P',t] if ca==0 else [t,'P']
                            elif t==middle:w+=['P','Q'] if cm==0 else ['Q','P']
                            elif t==right:w+=['Q',t] if cc==0 else [t,'Q']
                            else:w.append(t)
                        necks.append(w)
                    for s1,n1 in contract(necks,'P',group):
                        for s2,n2 in contract(n1,'Q',group):
                            if any(len(w)==1 for w in n2):continue
                            factor=N**sum(not w for w in n2);n2=[w for w in n2 if w]
                            written=[t for w in n2 for t in w]
                            base=Fraction(-c0*sf*(-1)**(ca+cm+cc)*factor*s1*s2*
                                _reorder_sign(written,ref,parity),6)
                            for c,dD,dE in jets:
                                LD,LE=dD+(SD,),dE+(SE,);cf=base*c;ws=[]
                                for w in n2:
                                    s,cw=canon_word(tuple(LD if t==left else LE if t==right else letters[t] for t in w))
                                    if not s:cf=0;break
                                    cf*=s;ws.append(cw)
                                if cf:
                                    s,v=canon_multi(ws)
                                    if s:out[v]+=s*cf
    return {w:c for w,c in out.items() if c}
