"""Replay corrected dMK (2.20) in free supercommutative graviton generator symbols."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys
from os import environ
from hashlib import sha256
vend = Path(environ['DMK_VENDOR']).resolve()
assert sha256((vend/'fortuity/bmn/dmk.py').read_bytes()).hexdigest() == '535ec58677878ffe5b613ec93bc4539c7e58eeefc8e532b98e791ccff0fc86b8'
assert sha256((vend/'fortuity/bmn/words.py').read_bytes()).hexdigest() == 'f74fba895060d3d5fb88922fff4ff38ff8d03546a8d56e7b25700fa19e93fe3e'
sys.path.insert(0,str(vend))
from fortuity.bmn import dmk

def atom(name,idx):
    return {((name,idx),):F(1)}

def to_monomial(w):
    odd=[a for a in w if a[0].startswith('v')]
    if len(odd)!=len(set(odd)):return (),0
    return tuple(sorted(w)),(-1)**sum(a>b for i,a in enumerate(odd) for b in odd[i+1:])

def mul(A,B):
    out=defaultdict(F)
    for a,b in product(A,B):
        k,s=to_monomial(a+b)
        if s: out[k]+=A[a]*B[b]*s
    return {k:c for k,c in out.items() if c}

def u2(i,j):
    return atom('u2',tuple(sorted((i,j))))
def u3(i,j,k):
    return atom('u3',tuple(sorted((i,j,k))))
def rawv2(i,j):
    return atom('v2',(i,j))
def v2(i,j):
    return dict(dmk.add((1,rawv2(i,j)), *[(F(-1,3),rawv2(a,a)) for a in range(3) if i==j]))
def rawv3(i,j,k):
    return atom('v3',(*sorted((i,j)),k))
def v3(i,j,k):
    return dict(dmk.add((1,rawv3(i,j,k)), *[(F(-1,4),rawv3(j,a,a)) for a in range(3) if i==k], *[(F(-1,4),rawv3(i,a,a)) for a in range(3) if j==k]))

dmk.mul_expr=mul
for name in ('u2','u3','v2','v3'):
    g=lru_cache(None)(globals()[name])
    setattr(dmk,name,g);dmk._RAW[name]=g
assert not dmk.relation_2_20((0,0,0))
total_separate=sum(len(dmk.term_expr(f,dict(i=0,j=0,k=0),use_R=True)) for _,_,f in dmk.TERMS)
assert total_separate > 0
coef,label,f=dmk.TERMS[0]
assert dmk.term_expr(f,dict(i=0,j=0,k=0),use_R=True)
print({'Eq_2_20_highest_component':'zero in free graded-commutative component polynomial ring', 'summands':len(dmk.TERMS),'separate_summand_monomials':total_separate,'negative_control':label})
