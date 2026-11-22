"""Exact transpose of the twice-iterated one-loop map into a pure gaugino power."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import argparse
import json
from pathlib import Path

import mixed_incoming as traces
from fortuity.full.letters import Q_letter, is_ferm
from fortuity.full.words import canon_multi
from orthogonal_glueball_first_pair import GHOSTS, PHI, orientations, remove_once, single_quantum
from orthogonal_odd_trace_primitive import primitive
from rank_four_lowered_replay import Q1d_multi, terminal


S = traces.so_canon((GHOSTS,))[1][0]


def tree_split(word):
    """Return a homotopy and the sorted section, with Q0 h(w)=w-section(w)."""
    if any(len(w)%2 for w in word):
        return primitive(word), {}
    work = [list(t) for t in word]
    image = Counter()
    permutation = 1
    def merge(k, i):
        x,y=work[k][i:i+2]
        jet=(x[0]+y[0], x[1]+y[1], tuple(sorted(x[2]+y[2])))
        coeff=sum(c for c,a,b in Q_letter(jet) if (a,b)==(x,y))
        assert coeff in (1,-1)
        before=sum(is_ferm(v) for tr in work[:k] for v in tr)+sum(map(is_ferm,work[k][:i]))
        raw=tuple(tuple(t) for t in work[:k]+[work[k][:i]+[jet]+work[k][i+2:]]+work[k+1:])
        s,w=traces.so_canon(raw)
        if s:
            image[w]+=permutation*((-1)**before)*s*coeff
    for k, trace in enumerate(work):
        for end in range(len(trace)-1,0,-1):
            for i in range(end):
                x,y=trace[i:i+2]
                if x>y:
                    merge(k,i)
                    trace[i],trace[i+1]=y,x
                    if is_ferm(x) and is_ferm(y):
                        permutation=-permutation
        for i,(x,y) in enumerate(zip(trace,trace[1:])):
            if x==y and is_ferm(x):
                merge(k,i)
                return {w:c for w,c in image.items() if c}, {}
    orientation,w=traces.so_canon(tuple(tuple(t) for t in work))
    return ({key:c for key,c in image.items() if c},
            {w:permutation*orientation} if orientation else {})


@lru_cache(None)
def first_score(colours, word):
    odd=[w for w in word if len(w)%2]
    if len(odd)!=1 or any(w!=S and w!=odd[0] for w in word):
        return 0
    active=odd[0]
    if len(active) not in (3,5):
        return 0
    n=(len(active)-1)//2
    sign,w=terminal((S,)*n)
    return sign*single_quantum(colours,active).get(w,0)


def predecessors(intermediates):
    support=set()
    def store(raw):
        s,w=canon_multi(raw)
        if s:
            support.add(w)
    templates={}
    for word in intermediates:
        for n in (1,2):
            for out in set(combinations(word,n)):
                letters=tuple(l for t in out for l in t)
                if not all(x in letters for x in GHOSTS):
                    continue
                initial=tuple(sorted(remove_once(letters,GHOSTS)+PHI))
                rest=remove_once(word,out)
                templates[rest,initial]=n==1 or templates.get((rest,initial),False)
    for (rest,initial), joined in templates.items():
        for trace in orientations(initial):
            store(rest+(trace,))
        if joined:
            for length in range(2,len(initial)//2+1):
                for first in set(combinations(initial,length)):
                    for left in orientations(first):
                        for right in orientations(remove_once(initial,first)):
                            store(rest+(left,right))
    return sorted(support)


def tree_columns(support):
    columns=set()
    for w in support:
        for k,tr in enumerate(w):
            for i in range(len(tr)):
                rot=tr[i:]+tr[:i]
                x,y=rot[:2]
                if set(x[2]) & set(y[2]):
                    continue
                jet=x[0]+y[0],x[1]+y[1],tuple(sorted(x[2]+y[2]))
                if not any((a,b)==(x,y) for _,a,b in Q_letter(jet)):
                    continue
                s,v=traces.so_canon(w[:k]+((jet,)+rot[2:],)+w[k+1:])
                if s:
                    columns.add(v)
    return columns


def quantum_columns(support):
    columns=set()
    for w in support:
        for size in (1,2):
            for subset in set(combinations(w,size)):
                leaves=tuple(x for t in subset for x in t)
                if any(g not in leaves for g in GHOSTS):
                    continue
                left=remove_once(leaves,GHOSTS)
                rest=remove_once(w,subset)
                for i in range(3):
                    psi=(0,0,tuple(j for j in range(3) if j!=i))
                    initial=left+(psi,PHI[i])
                    for trace in orientations(initial):
                        s,v=traces.so_canon(rest+(trace,))
                        if s:
                            columns.add(v)
                    if size==1:
                        for length in range(2,len(initial)//2+1):
                            for cut in set(combinations(initial,length)):
                                if (psi in cut)==(PHI[i] in cut) and PHI[i] not in left:
                                    continue
                                for a in orientations(cut):
                                    for b in orientations(remove_once(initial,cut)):
                                        s,v=traces.so_canon(rest+(a,b))
                                        if s:
                                            columns.add(v)
    return columns


def first_factorization(colours,power,first):
    columns=tree_columns(first) | quantum_columns([(S,)*power])
    sign,wanted=terminal((S,)*power)
    for col in columns:
        left=Fraction()
        for v,c in Q1d_multi(col,colours,'SO_odd').items():
            s,w=terminal(v)
            if s and w==wanted:
                left+=s*sign*c
        right=sum(c*first.get(w,0) for w,c in traces.q({col:1}).items())
        assert left==right,(col,left,right)
    return len(columns)


def second_factorization(colours,first,second):
    columns=tree_columns(second) | quantum_columns(first)
    for v in columns:
        left=Fraction()
        for raw,n in Q1d_multi(v,colours,'SO_odd').items():
            sig,w=traces.so_canon(raw)
            if sig:
                left-=n*sig*first.get(w,0)
        right=sum(c*second.get(w,0) for w,c in traces.q({v:1}).items())
        assert left==right,(v,left,right)
    return len(columns)


def second_dual(colours,power):
    intermediate=predecessors([(S,)*power])
    first={w:Fraction(first_score(colours,w)) for w in intermediate if first_score(colours,w)}
    upper=predecessors(first)
    scores=[]
    second={}
    for w in upper:
        h,section=tree_split(w)
        expected=Counter({w:1})
        expected.subtract(section)
        assert traces.q(h)=={v:n for v,n in expected.items() if n}
        value=Fraction()
        for pre,c in h.items():
            if not any(len(x[2])==2 for t in pre for x in t):
                continue
            for raw,n in Q1d_multi(pre,colours,'SO_odd').items():
                sig,v=traces.so_canon(raw)
                if sig and v in first:
                    value-=c*n*sig*first[v]
        if value:
            second[w]=value
            encoded=[[(PHI+GHOSTS).index(l) for l in tr] for tr in w]
            scores.append([encoded,str(value)])
    return dict(colours=colours,power=power,first_support=len(first),
                first_factorization_columns=first_factorization(colours,power,first),
                second_factorization_columns=second_factorization(colours,first,second),
                tested_words=len(upper),scores=scores)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',nargs=2,type=int,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not __debug__ or min(a.case)<3:
        p.error('Run with assertions enabled on a nontrivial colour algebra and power')
    d=second_dual(*a.case)
    a.out.write_text(json.dumps(d,separators=(',',':'))+'\n')
    print(d['first_support'],d['tested_words'],len(d['scores']))
