"""SO(8) trace graphs and single-epsilon graphs, including graded signs.

An invariant is (paths, traces). There are zero or four paths. Four paths
are contracted in their displayed order with epsilon_(1...8); coefficients
of paths precede coefficients of traces. This is the raw epsilon contraction,
so four identical skew even matrices give 2^4*4! times the Pfaffian.
Two epsilon tensors reduce to deltas; these graphs span both SO(8) parities.
"""
from functools import lru_cache
from collections import defaultdict
from itertools import product
from fortuity.bmn.words import canon_word, canon_multi, QLET, is_ferm, contents_of_sector
from fortuity.hwcoh.wordenum import multitraces_with_content


def ferm(w):
    return sum(is_ferm(l) for l in w)


@lru_cache(None)
def trace_canon(w):
    s,c=canon_word(w)
    if not s:return 0,None
    sr,cr=canon_word(w[::-1])
    sr *= (-1)**(len(w)+ferm(w)*(ferm(w)-1)//2)
    if c==cr and s!=sr:return 0,None
    return min(((c,s),(cr,sr)))[1],min(c,cr)


@lru_cache(None)
def path_canon(w):
    rev=w[::-1]
    sr=(-1)**(len(w)+1+ferm(w)*(ferm(w)-1)//2)
    if w==rev and sr==-1:return 0,None
    return (sr,rev) if rev<w else (1,w)


def canon(paths,traces):
    sign=1;ps=[];ts=[]
    for w in paths:
        s,c=path_canon(tuple(w))
        if not s:return 0,None
        sign*=s;ps.append(c)
    for w in traces:
        s,c=trace_canon(tuple(w))
        if not s:return 0,None
        sign*=s;ts.append(c)
    s,ps=canon_multi(ps)
    if not s:return 0,None
    sign*=s
    s,ts=canon_multi(ts)
    if not s:return 0,None
    return sign*s,(ps,ts)


def count(inv):
    c=[0]*7
    for w in inv[0]+inv[1]:
        for l in w:c[l]+=1
    return tuple(c)


@lru_cache(None)
def open_words(c):
    out=set();remaining=list(c);w=[]
    def rec():
        if len(w)==sum(c):
            s,cw=path_canon(tuple(w))
            if s:out.add(cw)
            return
        for l,n in enumerate(remaining):
            if n:
                remaining[l]-=1;w.append(l);rec();w.pop();remaining[l]+=1
    rec()
    return tuple(sorted(out))


@lru_cache(None)
def trace_monomials(c):
    if not any(c):return (((),()),)
    out=set()
    for m in multitraces_with_content(c):
        s,inv=canon((),m)
        if s:out.add(inv)
    return tuple(sorted(out))


@lru_cache(None)
def epsilon_paths(c):
    if sum(c)<4:return ()
    candidates=[]
    for sub in product(*(range(v+1) for v in c)):
        if 1<=sum(sub)<=sum(c)-3:
            for w in open_words(sub):candidates.append((w,sub,ferm(w)%2))
    candidates.sort();out=[]
    def rec(rem,start,chosen):
        if len(chosen)==4:
            if not any(rem):out.append(tuple(chosen))
            return
        if sum(rem)<4-len(chosen):return
        for i in range(start,len(candidates)):
            w,cc,odd=candidates[i]
            if all(a<=b for a,b in zip(cc,rem)):
                rec(tuple(b-a for a,b in zip(cc,rem)),i+odd,chosen+[w])
    rec(c,0,[])
    return tuple(out)


@lru_cache(None)
def invariants_content(c,parity):
    if parity==0:return trace_monomials(c)
    if sum(c)<4:return ()
    out=set()
    for pc in product(*(range(v+1) for v in c)):
        if sum(pc)<4:continue
        tc=tuple(a-b for a,b in zip(c,pc))
        if sum(tc)==1:continue
        traces=trace_monomials(tc)
        if not traces:continue
        for ps in epsilon_paths(pc):
            for _,ts in traces:out.add((ps,ts))
    return tuple(sorted(out))


def sector_basis(sec,parity):
    return [w for c in contents_of_sector(sec) for w in invariants_content(c,parity)]


@lru_cache(None)
def action(inv,which='Q'):
    out=defaultdict(int);par=0
    paths,traces=inv
    for typ,seq in enumerate((paths,traces)):
        for i,w in enumerate(seq):
            for j,l in enumerate(w):
                if which=='Q':reps=QLET.get(l,())
                else:
                    m=int(which);reps=[]
                    if l<3:
                        for k in range(3):
                            if len({m,l,k})==3:
                                e=1 if (m,l,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1
                                reps.append((e,(3+k,)))
                    elif l==3+m:reps=[(1,(6,))]
                for c,rep in reps:
                    new=seq[:i]+(w[:j]+rep+w[j+1:],)+seq[i+1:]
                    s,iv=canon(new,traces) if typ==0 else canon(paths,new)
                    if s:out[iv]+=(-1)**par*c*s
                par+=is_ferm(l)
    return {w:c for w,c in out.items() if c}


def action_expr(expr,which='Q'):
    out=defaultdict(int)
    for w,c in expr.items():
        for v,d in action(w,which).items():out[v]+=c*d
    return {w:c for w,c in out.items() if c}


def clear_caches():
    for fn in (open_words,epsilon_paths,invariants_content,trace_monomials,action):fn.cache_clear()
