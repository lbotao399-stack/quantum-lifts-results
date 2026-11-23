#!/usr/bin/env python3
"""Exact ordered-trace coefficients of new SO13 scalar relation candidates.

The generating function is Pf(J+A)/Pf(J), expanded as exp(1/2 tr log).
The auxiliary matrices are precisely upstream H,E,F. Differential marker
polynomials are written explicitly below, including their integer factors.
At total physical length 14 or 15 every coefficient vanishes in 13 colours.
One-fermion coefficients use one additional odd marker, hence no odd/odd
exchange is hidden in the trace expansion.
"""
import argparse,json,time,random
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations,product
from math import factorial,prod
from pathlib import Path

ZERO=(0,)*18
def const(c=1):return {ZERO:c}
def variable(k):
    a=list(ZERO);a[k]=1;return {tuple(a):1}
def add(*polys):
    out=Counter()
    for p in polys:out.update(p)
    return {m:c for m,c in out.items() if c}
def scale(p,c):return {m:v*c for m,v in p.items() if v*c}
def mul(p,q):
    out=Counter()
    for m,c in p.items():
        for n,d in q.items():out[tuple(a+b for a,b in zip(m,n))]+=c*d
    return {m:c for m,c in out.items() if c}
def power(p,k):
    out=const()
    for _ in range(k):out=mul(out,p)
    return out
def dot(a,b):
    # Inverse of (1/2)tr(M_a M_b), for M=(H,E,F).
    return add(mul(a[0],b[0]),scale(mul(a[1],b[2]),2),scale(mul(a[2],b[1]),2))

X=[[variable(3*i+j) for j in range(3)] for i in range(3)]
Y=[[variable(9+3*i+j) for j in range(3)] for i in range(3)]
D=add(*(scale(mul(mul(X[0][p[0]],X[1][p[1]]),X[2][p[2]]),(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))) for p in permutations(range(3))))
cross=[add(mul(X[0][(i+1)%3],X[1][(i+2)%3]),scale(mul(X[0][(i+2)%3],X[1][(i+1)%3]),-1)) for i in range(3)]
# Twice the inverse-metric image of the ordinary cross covector.
Z=[scale(cross[0],2),cross[2],cross[1]]
PAIR=add(*(dot(X[i],Y[i]) for i in range(3)))
AA=dot(X[0],X[0]);ZZ=dot(Z,Z)

def marker_polynomials(shape):
    if shape==(5,5,5):
        return [mul(power(D,4),PAIR)],power(D,5)
    assert shape==(7,5,3)
    sources=[mul(power(D,2),mul(PAIR,mul(AA,ZZ))),
             mul(power(D,2),mul(add(*(mul(dot(X[0],X[i]),dot(X[0],Y[i])) for i in range(3))),ZZ)),
             mul(power(D,2),mul(add(*(mul(dot(Z,X[i]),dot(Z,Y[i])) for i in range(3))),AA))]
    return sources,mul(power(D,3),mul(AA,ZZ))

def code(m):
    assert sum(m[9:])<=1
    return sum(n<<(4*i) for i,n in enumerate(m[:9])),next((i for i,n in enumerate(m[9:]) if n),-1)

@lru_cache(None)
def necklace(w):
    if len(w)<2:return 0,()
    a=min(w[i:]+w[:i] for i in range(len(w)))
    r=w[::-1];b=min(r[i:]+r[:i] for i in range(len(r)))
    if a==b and len(w)%2:return 0,()
    return ((-1)**len(w),b) if b<a else (1,a)

def canonical(ws):
    out=[];s=1
    for w in ws:
        t,u=necklace(tuple(w))
        if not t:return 0,()
        s*=t;out.append(u)
    return s,tuple(sorted(out))

ENTRIES=(((0,0,1),(1,1,-1)),((0,1,1),),((1,0,1),))

@lru_cache(None)
def auxtrace(w):
    states={(0,0,0,-1):1,(1,1,0,-1):1}
    for l in w:
        out=Counter()
        for (start,end,x,y),c in states.items():
            for a,entries in enumerate(ENTRIES):
                if l>=3:
                    if y>=0:continue
                    nx,ny=x,3*(l-3)+a
                else:nx,ny=x+(1<<(4*(3*l+a))),y
                for i,j,t in entries:
                    if end==i:out[start,j,nx,ny]+=c*t
        states={k:v for k,v in out.items() if v}
    out=Counter()
    for (start,end,x,y),c in states.items():
        if start==end:out[x,y]+=c
    return {k:v for k,v in out.items() if v}

@lru_cache(None)
def trace(w):
    assert necklace(w)==(1,w)
    orbit={v[i:]+v[:i] for v in (w,w[::-1]) for i in range(len(w))}
    factor=Fraction((-1)**(len(w)+1)*len(orbit),2*len(w))
    return {k:factor*v for k,v in auxtrace(w).items()}

def times(p,q):
    out=Counter()
    for (x,y),c in p.items():
        for (v,z),d in q.items():
            if y>=0 and z>=0:continue
            out[x+v,max(y,z)]+=c*d
    return {k:v for k,v in out.items() if v}

@lru_cache(maxsize=20000)
def traces_product(w):
    if not w:return {(0,-1):Fraction(1)}
    return times(traces_product(w[:-1]),trace(w[-1]))

class Coefficient:
    def __init__(self,p):
        self.raw=p
        self.weights={code(m):c*prod(factorial(n) for n in m) for m,c in p.items()}
    @lru_cache(None)
    def __call__(self,w):
        s,w=canonical(w)
        if not s:return Fraction()
        multiplicity=prod(factorial(n) for n in Counter(w).values())
        p=traces_product(w)
        return Fraction(s*sum(c*p.get(m,0) for m,c in self.weights.items()),multiplicity)

PAIRS=((1,2),(2,0),(0,1))
def differential(w):
    out=Counter()
    for i,t in enumerate(w):
        for j,l in enumerate(t):
            if l<3:continue
            a,b=PAIRS[l-3]
            for pair,c in [((a,b),1),((b,a),-1)]:
                s,v=canonical(w[:i]+(t[:j]+pair+t[j+1:],)+w[i+1:])
                if s:out[v]+=s*c
    return {k:v for k,v in out.items() if v}

@lru_cache(None)
def predecessors(w):
    candidates=set()
    for i,t in enumerate(w):
        for j in range(len(t)):
            r=t[j:]+t[:j];a,b=r[:2]
            if a==b:continue
            l=3+next(k for k,pair in enumerate(PAIRS) if {a,b}==set(pair))
            s,v=canonical(w[:i]+((l,)+r[2:],)+w[i+1:])
            if s:candidates.add(v)
    return {v:differential(v)[w] for v in candidates if w in differential(v)}

def image(source,w):return sum(c*source(v) for v,c in predecessors(w).items())
def decode(s):return tuple(tuple('ABCabc'.index(l) for l in t) for t in s.split('.'))
def encode(w):return '.'.join(''.join('ABCabc'[l] for l in t) for t in w)

def matrix_rank(rows):
    rows=[list(map(Fraction,r)) for r in rows];rank=0
    if not rows:return 0
    for j in range(len(rows[0])):
        pivot=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        q=rows[rank][j];rows[rank]=[x/q for x in rows[rank]]
        for i in range(rank+1,len(rows)):
            q=rows[i][j]
            if q:rows[i]=[x-q*y for x,y in zip(rows[i],rows[rank])]
        rank+=1
    return rank

def solve(a,b):
    n=len(a);rows=[list(map(Fraction,r))+[Fraction(v)] for r,v in zip(a,b)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[pivot]=rows[pivot],rows[j]
        c=rows[j][j];rows[j]=[x/c for x in rows[j]]
        for i in range(n):
            if i==j:continue
            c=rows[i][j]
            if c:rows[i]=[x-c*y for x,y in zip(rows[i],rows[j])]
    return [r[-1] for r in rows]

def determinant(rows):
    rows=[list(map(Fraction,r)) for r in rows];det=Fraction(1)
    for j in range(len(rows)):
        pivot=next(i for i in range(j,len(rows)) if rows[i][j])
        if pivot!=j:rows[j],rows[pivot]=rows[pivot],rows[j];det=-det
        c=rows[j][j];det*=c
        for i in range(j+1,len(rows)):
            q=rows[i][j]/c
            rows[i]=[x-q*y for x,y in zip(rows[i],rows[j])]
    return det

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--shape',type=int,nargs=3,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic();shape=tuple(args.shape)
    raw,terminal=marker_polynomials(shape);sources=list(map(Coefficient,raw));terminal=Coefficient(terminal)
    probes=['AA.AA.ABC.BB.BB.CC.CC','AA.ABC.ABC.ABC.BB.CC','ABC.ABC.ABC.ABC.ABC'] if shape==(5,5,5) else [
        'AA.AA.AA.ABC.BB.BB.CC','AA.AA.ABC.ABC.ABC.BB','AA.ABC.ABC.ABC.ABCB']
    # The final two probes for (7,5,3) are generated deterministically below.
    rng=random.Random(20261123+sum((i+1)*n for i,n in enumerate(shape)))
    layouts=[(2,2,2,2,2,2,3),(2,2,2,2,3,4),(2,2,2,3,3,3),(2,2,2,2,2,5),
             (2,2,3,3,5),(2,2,2,4,5),(3,3,3,3,3),(2,2,2,2,7),(2,2,3,4,4)]
    samples=[i for i,n in enumerate(shape) for _ in range(n)]
    for _ in range(120):
        rng.shuffle(samples);layout=rng.choice(layouts);w=[];j=0
        for n in layout:w.append(tuple(samples[j:j+n]));j+=n
        s,w=canonical(w)
        if s:probes.append(encode(w))
    rows=[];selected=[];selected_values=[];seen=set()
    for probe in probes:
        s,w=canonical(decode(probe))
        if not s or w in seen or tuple(sum(t.count(i) for t in w) for i in range(3))!=shape:continue
        seen.add(w)
        values=[image(source,w) for source in sources]+[terminal(w)]
        row={'probe':encode(w),'Q0_source_coefficients':list(map(str,values[:-1])),'terminal_coefficient':str(values[-1])}
        row['basis_probe']=matrix_rank(selected_values+[values])>len(selected_values)
        rows.append(row)
        if row['basis_probe']:
            selected_values.append(values);selected.append(row)
            print(json.dumps(row),flush=True)
        if len(selected_values)==len(sources)+1:break
    result={'shape':shape,'colours':13,'field':'Q','source_marker_terms':[len(p) for p in raw],
            'terminal_marker_terms':len(terminal.raw),'probes':rows,'seconds':time.monotonic()-start,
            'scope':'actual finite-colour Pfaffian identities; probe values only, not yet a quotient or novelty certificate'}
    if len(selected_values)==len(sources)+1:
        dual=solve(list(zip(*selected_values)),[0]*len(sources)+[1])
        result.update(basis_probe_words=[r['probe'] for r in selected],
                      exact_probe_matrix=[[str(c) for c in r] for r in selected_values],
                      exact_determinant=str(determinant(selected_values)),
                      dual_coefficients=list(map(str,dual)),
                      dual_pairing_on_sources=[str(sum(dual[i]*selected_values[i][j] for i in range(len(dual)))) for j in range(len(sources))],
                      dual_pairing_on_terminal=str(sum(dual[i]*selected_values[i][-1] for i in range(len(dual)))))
        print('FULL_RANK',shape,'det',result['exact_determinant'],'dual',result['dual_coefficients'],flush=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
