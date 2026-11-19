"""Exact unprojected three-colour singlet dimensions for fixed letter multiplicities.

Use simple-root coordinates of A2. The adjoint has six roots and two zero
weights. Even and odd species contribute symmetric and exterior powers,
respectively. A Weyl denominator extracts the trivial isotypic component.
An independent constant-term check uses the two-sided denominator, divided
by the Weyl group order. Charge conjugation is not imposed.
"""
import argparse
from collections import Counter
from functools import lru_cache
from pathlib import Path
import replay_formal_blocks as s
from fortuity.full.letters import is_ferm

POS=((1,0),(0,1),(1,1));ROOTS=POS+tuple((-i,-j) for i,j in POS);ADJ=ROOTS+((0,0),(0,0))

def mul(a,b):
    result=Counter()
    for (i,j),n in a.items():
        for (k,l),m in b.items():result[i+k,j+l]+=n*m
    return {u:n for u,n in result.items() if n}

def denominator(roots):
    result={(0,0):1}
    for i,j in roots:result=mul(result,{(0,0):1,(i,j):-1})
    return result

@lru_cache(None)
def power(n,odd):
    values=[Counter({(0,0):1})]+[Counter() for i in range(n)]
    for i,j in ADJ:
        for d in (range(n,0,-1) if odd else range(1,n+1)):
            for (k,l),v in values[d-1].items():values[d][i+k,j+l]+=v
    return values[-1]

@lru_cache(None)
def size(profile):
    char={(0,0):1}
    for n,odd in profile:char=mul(char,power(n,odd))
    def ct(roots):return sum(n*char.get((-i,-j),0) for (i,j),n in denominator(roots).items())
    n=ct(POS);assert n>=0 and ct(ROOTS)==6*n;return n

def dimension(content):
    return size(tuple(sorted((v,bool(is_ferm(s.ALPH[i]))) for i,v in Counter(content).items())))
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('contents',nargs='+');b=a.parse_args();print({x:dimension(bytes.fromhex(x)) for x in b.contents})
