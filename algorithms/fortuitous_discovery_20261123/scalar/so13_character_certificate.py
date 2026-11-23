#!/usr/bin/env python3
"""Self-contained exact SO13 relation characters, independent of runtime paths.

Orthogonal invariant tensors have even-row diagrams. At degrees 14 and15
the omitted defining-vector diagrams are (2^14), and (2^15),(4,2^13).
Transposition turns these into (14,14), and (15,15),(14,14,1,1).
The transpose sign cancels the antisymmetric-adjoint projector sign, giving
the positive plethystic sum implemented in relation_trace().
"""
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product,combinations
from math import comb,factorial,prod

@lru_cache(None)
def partitions(n,minimum=1):
    if n==0:return ((),)
    return tuple((k,)+tail for k in range(minimum,n+1) for tail in partitions(n-k,k))

def transpose(shape):return tuple(sum(x>=j for x in shape) for j in range(1,max(shape,default=0)+1))

@lru_cache(None)
def representation_dimension(shape):
    conjugate=transpose(shape)
    return factorial(sum(shape))//prod(shape[i]-j+conjugate[j]-i-1 for i in range(len(shape)) for j in range(shape[i]))

@lru_cache(None)
def character(shape,cycles):
    if not cycles or cycles[-1]==1:return representation_dimension(shape)
    n=len(shape);length=cycles[-1];beads=tuple(x+n-i-1 for i,x in enumerate(shape));value=0
    for i,b in enumerate(beads):
        destination=b-length
        if destination<0 or destination in beads:continue
        sign=(-1)**sum(destination<v<b for v in beads)
        moved=sorted(beads[:i]+beads[i+1:]+(destination,),reverse=True)
        smaller=tuple(b-n+j+1 for j,b in enumerate(moved) if b-n+j+1)
        value+=sign*character(smaller,cycles[:-1])
    return value

def centralizer(cycles):return prod(k**m*factorial(m) for k,m in Counter(cycles).items())

@lru_cache(None)
def relation_trace(cycles):
    degree=sum(cycles)
    if degree<=13:return 0
    shapes={14:((14,14),),15:((15,15),(14,14,1,1))}[degree]
    types=sorted(Counter(cycles).items());value=0
    for crossed in product(*(range(m+1) for _,m in types)):
        vector_cycles=[];coefficient=1
        for (k,m),a in zip(types,crossed):
            vector_cycles += [k]*(2*(m-a))+[2*k]*a
            coefficient*=comb(m,a)
        vector_cycles=tuple(sorted(vector_cycles))
        value+=coefficient*sum(character(shape,vector_cycles) for shape in shapes)
    result,remainder=divmod(value,2**len(cycles));assert not remainder
    return result

def flavour(shape,distinguished=0):
    value=sum(Fraction(character(shape,c)*relation_trace(tuple(sorted(c+(1,)*distinguished))),centralizer(c)) for c in partitions(sum(shape)))
    assert value.denominator==1 and value>=0
    return int(value)

def predecessor_shapes(shape):
    out=[]
    for removed in combinations(range(3),2):
        s=tuple(n-(i in removed) for i,n in enumerate(shape))
        if min(s)>=0 and list(s)==sorted(s,reverse=True):out.append(tuple(n for n in s if n))
    return out

def certificate():
    rows=[]
    for shape in [(5,5,5),(7,5,3)]:
        source={','.join(map(str,s)):flavour(s,1) for s in predecessor_shapes(shape)}
        target=flavour(shape)
        rows.append(dict(shape=shape,source_Pieri_multiplicities=source,source_dimension=sum(source.values()),terminal_dimension=target))
    assert [(r['source_dimension'],r['terminal_dimension']) for r in rows]==[(1,2),(3,4)]
    return dict(passed=True,arithmetic='Q',omitted_defining_vector_diagrams={'14':[[2]*14],'15':[[2]*15,[4]+[2]*13]},
                transposed_diagrams={'14':[[14,14]],'15':[[15,15],[14,14,1,1]]},records=rows)

if __name__=='__main__':print(json.dumps(certificate(),indent=2))
