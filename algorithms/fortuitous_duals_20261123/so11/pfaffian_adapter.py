#!/usr/bin/env python3
"""Load the pinned Pfaffian coefficient class without unrelated census imports.

The coefficient class, constants, and scalar_code function are compiled
unchanged from the pinned upstream AST. Only graded SO trace canonicalization
is supplied locally, with the same tuple-letter order. No physical coefficient
or factorial is changed.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations,product
from math import factorial,prod
from pathlib import Path
from types import SimpleNamespace
import ast,hashlib,os

HERE=Path(__file__).resolve().parent
SOURCE=Path(os.environ.get('SO11_PFAFFIAN_PATH',HERE/'inputs/characters/orthogonal_glueball_first_pair.py'))

def odd(letter):
    return (letter[0]+letter[1])%2

@lru_cache(None)
def canonical_trace(word):
    if len(word)<2:return 0,()
    count=sum(map(odd,word));variants={}
    for source,initial in [(word,1),(word[::-1],(-1)**(len(word)+count*(count-1)//2))]:
        for j in range(len(word)):
            cut=sum(map(odd,source[:j]));orientation=initial*(-1)**(cut*(count-cut))
            name=source[j:]+source[:j]
            variants.setdefault(name,set()).add(orientation)
    name=min(variants)
    if len(variants[name])>1:return 0,()
    return next(iter(variants[name])),name

@lru_cache(None)
def canonical(words):
    out=[];orientation=1
    for word in words:
        local,name=canonical_trace(word)
        if not local:return 0,()
        orientation*=local;out.append(name)
    oddwords=[w for w in out if sum(map(odd,w))%2]
    if len(set(oddwords))<len(oddwords):return 0,()
    orientation*=(-1)**sum(x>y for i,x in enumerate(oddwords) for y in oddwords[i+1:])
    return orientation,tuple(sorted(out))

def sign(perm):
    return (-1)**sum(x>y for i,x in enumerate(perm) for y in perm[i+1:])

def _load():
    code=SOURCE.read_bytes();tree=ast.parse(code,filename=str(SOURCE))
    allowed={'PHI','PLUS','MINUS','GHOSTS','ENTRIES'}
    selected=[]
    for node in tree.body:
        if isinstance(node,(ast.FunctionDef,ast.ClassDef)) and node.name in {'scalar_code','PfaffianCoefficient'}:
            selected.append(node)
        elif isinstance(node,ast.Assign):
            assigned={x.id for t in node.targets for x in ast.walk(t) if isinstance(x,ast.Name)}
            if assigned and assigned<=allowed:selected.append(node)
    assert len(selected)==6,[(type(n).__name__,getattr(n,'name',None)) for n in selected]
    ns=dict(Counter=Counter,Fraction=Fraction,lru_cache=lru_cache,permutations=permutations,
            product=product,factorial=factorial,prod=prod,sign=sign,traces=SimpleNamespace(so_canon=canonical))
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(SOURCE),'exec'),ns)
    return ns,hashlib.sha256(code).hexdigest()

_namespace,SOURCE_SHA256=_load()
PfaffianCoefficient=_namespace['PfaffianCoefficient']
PHI=_namespace['PHI'];GHOSTS=_namespace['GHOSTS']

class ScalarPfaffian(PfaffianCoefficient):
    """The scalar det(d/dX)^5 coefficient of the same Pfaffian series."""
    def __init__(self):super().__init__(5,0)
    def coefficient(self,words):
        multiplicity=prod(factorial(n) for n in Counter(words).values())
        poly=self.traces_product(words)
        return Fraction(sum(c*poly.get((0,x,0),0) for x,c in self.weights.items()))/multiplicity

def encode(word):
    return tuple(tuple(dict(zip('ABC+-',PHI+GHOSTS))[c] for c in w) for w in word)

def oracle(pf,word):
    orientation,canonical_words=canonical(encode(word))
    return orientation*pf.coefficient(canonical_words) if orientation else Fraction()
