#!/usr/bin/env python3
"""Exact rational finite-SU(3) classical candidate; Python standard library only."""
import gzip
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parent
def clean(a):
    return {w:c for w,c in a.items() if c}
def odd(l):
    return (len(l[2])+1)%2
@lru_cache(None)
def trace_parity(w):
    return sum(map(odd,w))%2
@lru_cache(maxsize=200000)
def cyclic(w):
    if len(w)==1:
        return 0,()
    candidates={}
    s=1
    for _ in w:
        if w in candidates and candidates[w]!=s:
            return 0,()
        candidates[w]=s
        s*=(-1)**(odd(w[0])*trace_parity(w[1:]))
        w=w[1:]+w[:1]
    best=min(candidates)
    return candidates[best],best
def ordered(w,parity):
    exterior=[x for x in w if parity(x)]
    if len(set(exterior))!=len(exterior):
        return 0,()
    s=(-1)**sum(a>b for i,a in enumerate(exterior) for b in exterior[i+1:])
    return s,tuple(sorted(w))
def canonical(w):
    traces=[]
    s=1
    for tr in w:
        t,v=cyclic(tr)
        if not t:
            return 0,()
        s*=t
        traces.append(v)
    t,v=ordered(traces,trace_parity)
    return s*t,v
def add_term(a,w,c):
    s,v=canonical(w)
    if s:
        a[v]+=c*s
def multiply(a,b):
    r=Counter()
    for w,c in a.items():
        for v,d in b.items():
            add_term(r,w+v,c*d)
    return clean(r)
@lru_cache(None)
def q_letter(l):
    a,b,i=l
    out=[]
    for u,v in product(range(a+1),range(b+1)):
        for n in range(len(i)+1):
            for j in combinations(i,n):
                k=tuple(x for x in i if x not in j)
                if (u or v or j) and (a-u or b-v or k):
                    s=(-1)**(len(i)+(len(j)+1)*len(k)+sum(x>y for x in j for y in k))
                    out.append((s,((u,v,j),(a-u,b-v,k))))
    return tuple(out)
def operation(l,op):
    a,b,i=l
    x,y,z=op
    d=a+b+len(i)
    if x=='Q':
        return q_letter(l)
    if x=='D':
        if y<2:
            r=((a,b)[y]+1,(a+(y==0),b+(y==1),i))
        elif y-2 not in i:
            r=((-1)**sum(k>y-2 for k in i),(a,b,tuple(sorted(i+(y-2,)))))
        else:
            return ()
    elif x=='E' and z in i and y not in i:
        r=((-1)**sum(min(y,z)<k<max(y,z) for k in i),(a,b,tuple(sorted((set(i)-{z})|{y}))))
    elif x=='J' and b:
        r=(a+1,(a+1,b-1,i))
    elif x=='T' and (a,b)[y] and z not in i:
        r=((-1)**(len(i)+sum(k<z for k in i)),(a-(y==0),b-(y==1),tuple(sorted(i+(z,)))))
    elif x=='S' and y in i and d>1:
        r=((d-1)*(-1)**(len(i)+sum(k<y for k in i)),(a,b,tuple(k for k in i if k!=y)))
    elif x=='K' and (a,b)[y] and d>1:
        r=(d-1,(a-(y==0),b-(y==1),i))
    else:
        return ()
    return ((r[0],(r[1],)),)
def act(a,op,quotient=True):
    out=Counter()
    odd_op=op[0] in ('Q','S','T') or op[0]=='D' and op[1]>=2
    for w,c in a.items():
        prefix=0
        for k,tr in enumerate(w):
            for j,l in enumerate(tr):
                for d,repl in operation(l,op):
                    v=w[:k]+(tr[:j]+repl+tr[j+1:],)+w[k+1:]
                    t=c*d*(-1)**(prefix*odd_op)
                    if quotient:
                        add_term(out,v,t)
                    else:
                        out[v]+=t
                prefix+=odd(l)
    return clean(out)
@lru_cache(None)
def factor(g):
    plus,minus=(1,0,()),(0,1,())
    inside=(plus,)*g[0]+(minus,)*g[1]
    inside+=tuple((0,0,(i,)) for i,n in enumerate(g[2:5]) for _ in range(n))
    words=set(permutations(inside))
    a={(w,):F(-1 if plus in w and minus in w and w.index(minus)<w.index(plus) else 1,len(words)) for w in words}
    axes=[i+2 for i in range(3) if g[7]>>i&1]+[0]*g[5]+[1]*g[6]
    for i in axes:
        a=act(a,('D',i,0),quotient=False)
    b=Counter()
    for w,c in a.items():
        add_term(b,w,c)
    return clean(b)
def diagonality(l):
    return (l[2].bit_count()+1)%2
def slot_parity(w):
    return sum(map(diagonality,w))%2
@lru_cache(None)
def diag(tr):
    jets=tuple((a,b,sum(1<<i for i in th)) for a,b,th in tr)
    r=Counter()
    for colors in product(range(2),repeat=len(jets)):
        mult=int(not any(colors))+int(all(colors))+(-1)**len(jets)
        if not mult:
            continue
        left=tuple(x for x,k in zip(jets,colors) if not k)
        right=tuple(x for x,k in zip(jets,colors) if k)
        s,a=ordered(left,diagonality)
        t,b=ordered(right,diagonality)
        cross=sum(diagonality(x)*diagonality(y) for i,x in enumerate(jets) if colors[i] for j,y in enumerate(jets[i+1:],i+1) if not colors[j])
        if s*t:
            r[a,b]+=s*t*mult*(-1)**cross
    return clean(r)
def rho(a):
    out=Counter()
    for w,c in a.items():
        r={((),()):c}
        for tr in w:
            nxt=Counter()
            for (a,b),d in r.items():
                for (u,v),e in diag(tr).items():
                    s,x=ordered(a+u,diagonality)
                    t,y=ordered(b+v,diagonality)
                    if s*t:
                        nxt[x,y]+=d*e*s*t*(-1)**(slot_parity(b)*slot_parity(u))
            r=clean(nxt)
        out.update(r)
    return clean(out)

def main():
    archive=json.loads(gzip.decompress((ROOT/'su30_0.json.gz').read_bytes()))['selected_recipes']
    recipe=json.loads((ROOT/'R007_recipe.json').read_text())['source_linear_combination']
    factors={tuple(g) for r in archive for g in r}
    assert len(factors)==84 and all(not act(factor(g),('Q',0,0)) for g in factors)
    a=Counter()
    for i,c in recipe:
        b={():F(c)}
        for g in archive[i]:
            b=multiply(b,factor(tuple(g)))
        a.update(b)
    a=clean(a)
    sizes=[len(a)]
    for x,y in [(0,1),(1,1),(0,2),(1,2)]:
        a=act(a,('T',x,y))
        sizes.append(len(a))
    assert sizes==[348,1059,1864,4561,5080]
    charge=set()
    conjugate=Counter()
    for multi,c in a.items():
        letters=[l for tr in multi for l in tr]
        charge.add((sum(l[0] for l in letters),sum(l[1] for l in letters),
                    *(sum(i in l[2] for l in letters) for i in range(3)),
                    len(letters),sum(map(odd,letters))%2))
        exponent=sum(len(tr)+sum(map(odd,tr))*(sum(map(odd,tr))-1)//2 for tr in multi)
        add_term(conjugate,tuple(tr[::-1] for tr in multi),c*(-1)**exponent)
    assert charge=={(2,2,5,2,2,9,0)}
    assert clean(conjugate)=={w:-c for w,c in a.items()}
    ops=[('E',0,1),('E',1,2),('J',0,0)]
    ops += [('T',x,y) for x,y in product(range(2),range(3))]
    ops += [('S',i,0) for i in range(3)]+[('K',i,0) for i in range(2)]
    assert len(ops)==14
    for op in [('Q',0,0)]+ops:
        assert not act(a,op),op
    # Negative control: a raw archived product cannot replace the rational recipe.
    b={():1}
    for g in archive[recipe[0][0]]:
        b=multiply(b,factor(tuple(g)))
    for x,y in [(0,1),(1,1),(0,2),(1,2)]:
        b=act(b,('T',x,y))
    assert b and any(act(b,op) for op in ops)
    d=rho(a)
    canonical_data=[[w,str(c)] for w,c in sorted(d.items())]
    digest=hashlib.sha256(json.dumps(canonical_data,separators=(',',':')).encode()).hexdigest()
    assert digest=='89575f063dd2a39d7f45f3d7cd7a93dc1d2284141361aba79c29767e3b4e8ccd'
    assert len(d)==34104 and d[min(d)]==F(432,5)
    print('PASS: 5080 trace terms, Q0 and all 14 equations literal zero, 34104 Cartan coefficients match, first 432/5.')
if __name__=='__main__':
    main()
