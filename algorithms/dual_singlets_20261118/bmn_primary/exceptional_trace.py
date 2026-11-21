"""Adjoint permutation traces from finite exceptional root systems, over the integers."""
from collections import Counter
from functools import lru_cache


@lru_cache(None)
def root_data(group):
    rank=int(group[1:])
    if group not in ("E6","E7","E8"):
        raise ValueError("Expected a simply laced exceptional root system")
    cartan=[[int(i==j)*2 for j in range(rank)] for i in range(rank)]
    edges=[(0,1),(0,2),(2,3),(0,4)]+[(i,i+1) for i in range(4,rank-1)]
    for i,j in edges:
        cartan[i][j]=cartan[j][i]=-1
    simple=tuple(tuple(row) for row in cartan)
    known=set(simple)
    pending=list(simple)
    while pending:
        w=pending.pop()
        for i,root in enumerate(simple):
            if w[i]:
                v=tuple(a-w[i]*b for a,b in zip(w,root))
                if v not in known:
                    known.add(v)
                    pending.append(v)
    assert len(known)=={"E6":72,"E7":126,"E8":240}[group]
    return simple,tuple(sorted(known))


@lru_cache(maxsize=150000)
def chamber(group, shifted):
    simple,_=root_data(group)
    sign=1
    while 0 not in shifted:
        if min(shifted)>0:
            return sign,tuple(i-1 for i in shifted)
        k=next(i for i,n in enumerate(shifted) if n<0)
        shifted=tuple(a-shifted[k]*b for a,b in zip(shifted,simple[k]))
        sign=-sign
    return 0,()


@lru_cache(maxsize=40000)
def multiply(group, highest, k):
    simple,roots=root_data(group)
    result=Counter({highest:len(simple)})
    for root in roots:
        sign,other=chamber(group,tuple(a+1+k*b for a,b in zip(highest,root)))
        if sign:
            result[other]+=sign
    return tuple((w,n) for w,n in result.items() if n)


@lru_cache(None)
def decomposition(group, cycles):
    if not cycles:
        return {(0,)*len(root_data(group)[0]):1}
    result=Counter()
    for w,n in decomposition(group,cycles[:-1]).items():
        for other,m in multiply(group,w,cycles[-1]):
            result[other]+=m*n
    return {w:n for w,n in result.items() if n}


@lru_cache(None)
def invariant_trace(group,cycles):
    # Both products of Adams characters are self-dual, including in type E6.
    # Their character inner product picks out the trivial representation.
    a,b=[],[]
    left=right=0
    for k in sorted(cycles,reverse=True):
        if left<=right:
            a.append(k)
            left+=k
        else:
            b.append(k)
            right+=k
    first=decomposition(group,tuple(sorted(a)))
    second=decomposition(group,tuple(sorted(b)))
    if len(first)>len(second):
        first,second=second,first
    return sum(n*second.get(w,0) for w,n in first.items())


@lru_cache(None)
def folded_roots():
    from f4_singlet_characters import CARTAN,reflect,roots
    simple=tuple(enumerate(CARTAN["F4"]))
    # The last simple root is short in the convention with row alpha_i.
    short={simple[-1][1]}
    pending=list(short)
    while pending:
        w=pending.pop()
        for s in simple:
            v=reflect(w,s)
            if v not in short:
                short.add(v)
                pending.append(v)
    long=set(roots("F4"))-short
    assert len(short)==len(long)==24
    return tuple(sorted(short)),tuple(sorted(long))


@lru_cache(None)
def folded_step(highest,k):
    from f4_singlet_characters import chamber as folded_chamber
    short,long=folded_roots()
    result=Counter({highest:2 if k%2 else 6})
    weights=[(short,k,1)]
    if not k%2:
        weights.append((long,k//2,2))
    for roots,exponent,multiplicity in weights:
        for root in roots:
            sign,target=folded_chamber("F4",tuple(a+1+exponent*b for a,b in zip(highest,root)))
            if sign:
                result[target]+=sign*multiplicity
    return tuple((w,n) for w,n in result.items() if n)


@lru_cache(None)
def folded_decomposition(cycles):
    if not cycles:
        return {(0,0,0,0):1}
    result=Counter()
    for w,n in folded_decomposition(cycles[:-1]).items():
        for other,m in folded_step(w,cycles[-1]):
            result[other]+=n*m
    return {w:n for w,n in result.items() if n}


@lru_cache(None)
def e6_traces(cycles):
    direct=invariant_trace("E6",cycles)
    ordered=sorted(cycles,reverse=True)
    first=folded_decomposition(tuple(sorted(ordered[::2])))
    second=folded_decomposition(tuple(sorted(ordered[1::2])))
    twisted=sum(n*second.get(w,0) for w,n in first.items())
    if (direct-twisted)%2:
        raise ArithmeticError("Nonintegral outer-involution components")
    return (direct+twisted)//2,(direct-twisted)//2
