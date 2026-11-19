"""Exact character bounds after quotienting terminal traces by the classical supercommutators."""
import argparse
import json
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path


def add(x,y): return tuple(i+j for i,j in zip(x,y))

def sub(x,y): return tuple(i-j for i,j in zip(x,y))

def parity(g): return (g[0]+g[1])%2

@lru_cache(None)
def generators(w):
    if min(w)<0:return ()
    return tuple(sorted((g for g in product(*(range(min(v,1)+1) if i<2 else range(v+1) for i,v in enumerate(w))) if sum(g)>0 and sum(g)%2==0),key=lambda g:(sum(g),g)))

@lru_cache(None)
def words(w):
    if min(w)<0:return ()
    gen=generators(w)
    def each(rem,start):
        if not any(rem):yield ();return
        for i in range(start,len(gen)):
            g=gen[i];r=sub(rem,g)
            if min(r)<0:continue
            for t in each(r,i+parity(g)):yield (g,)+t
    return tuple(each(w,0))

def sign(perm):return (-1)**sum(x>y for i,x in enumerate(perm) for y in perm[i+1:])

def weight_terms(h):
    w=(h,h,4-h,4-h,4-h);rho=(2,1,0)
    for d,s in [(0,1),(1,-1)]:
        for perm in permutations(range(3)):
            yield (w[0]+d,w[1]-d)+tuple(w[2+i]+rho[i]-rho[perm[i]] for i in range(3)), s*sign(perm)

def audit_rules():
    from charge_maps import RUNTIME
    import sp_environment
    from fortuity.full.letters import Q_letter,is_ferm
    letters=[(1,0,()),(0,1,())]+[(0,0,(i,)) for i in range(3)]
    checks=[]
    for i,x in enumerate(letters):
        assert not Q_letter(x)
        for j in range(i,5):
            y=letters[j]
            if x==y and not is_ferm(x):continue
            e=(x[0]+y[0],x[1]+y[1],tuple(sorted(x[2]+y[2])));a={(l,r):c for c,l,r in Q_letter(e)}
            if x==y: assert a=={(x,y):1}
            else: assert len(a)==2 and abs(a[x,y])==1 and a[x,y]==-(-1)**(is_ferm(x)*is_ferm(y))*a[y,x]
            checks.append([i,j])
    assert len(checks)==12
    return checks

def compute():
    checks=audit_rules();result=[]
    for h in range(5):
        parts=[]
        for w,s in weight_terms(h):
            m=words(w);parts.append(dict(weight=w,sign=s,formal=len(m),cartan=sum(len(t)<=3 for t in m)))
        ub=sum(p['sign']*p['formal'] for p in parts);lb=sum(p['sign']*p['cartan'] for p in parts)
        assert 0<=lb<=ub
        result.append(dict(bidegree=[h,h],length=12-h,cartan=lb,formal_trace_quotient=ub,character_terms=parts))
    return dict(classical_letter_pairs=checks,rank=3,groups=['SO(7)','USp(6)'],blocks=result)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path);a=ap.parse_args();r=compute();
    if a.output:a.output.write_text(json.dumps(r,indent=2)+'\n')
    for t in r['blocks']:print(t['bidegree'],t['cartan'],t['formal_trace_quotient'])
