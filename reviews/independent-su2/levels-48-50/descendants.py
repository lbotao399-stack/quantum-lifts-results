"""Compare named BMN polynomials with the complete rational boundary span.

The variables and the three traces in O_n are those of KKLP (2.15).
Only the internal divided-jet coordinate for f differs, as specified in core().
"""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import permutations, product
from math import gcd, lcm
from pathlib import Path
import json
import sys

if not (Path(__file__).parent / 'su2_bmn_frame.py').exists():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import su2_bmn_frame as frame


def clean(row):
    return {k: v for k, v in row.items() if v}


def add_to(row, poly, scale=1):
    for key, value in poly.items():
        row[key] += scale * value


def multiply(left, right):
    result = defaultdict(int)
    for (a, e), u in left.items():
        for (b, f), v in right.items():
            sign = frame.wedge(e, f)
            if sign:
                result[a+b, e|f] += sign*u*v
    return clean(result)


def mult(*factors):
    value = {(0, 0): 1}
    for factor in factors:
        value = multiply(value, factor)
    return value


def eps(i, j, k):
    return ((i-j)*(j-k)*(k-i))//2


axes = range(3)
perms = tuple(permutations(axes))
# Four bits store the exponent of each of the twelve commuting components.
phi = tuple(tuple({(1 << (4*(3*i+j)), 0): 1} for j in axes) for i in axes)
psi = tuple(tuple({(0, 1 << (3*i+j)): 1} for j in axes) for i in axes)
field = tuple({(1 << (4*(9+j)), 0): 1} for j in axes)


def pair(x, y):
    result = defaultdict(int)
    for i,j,n in ((0,2,1),(1,1,2),(2,0,1)):
        add_to(result, multiply(x[i], y[j]), n)
    return clean(result)


def triple(x, y, z):
    result = defaultdict(int)
    for i,j,k in perms:
        add_to(result, mult(x[i],y[j],z[k]), -eps(i,j,k))
    return clean(result)


scalar_fermion = tuple(tuple(pair(phi[i], psi[j]) for j in axes) for i in axes)
field_fermion = tuple(pair(field, psi[i]) for i in axes)
fermion_pair = tuple(tuple(pair(x,y) for y in psi) for x in psi)
fermion_triple = tuple(tuple(tuple(triple(x,y,z) for z in psi) for y in psi) for x in psi)


@cache
def core(n):
    # KKLP (2.3) has Q_0 psi = -2i [phi,phi] cyclically and
    # Q_0 f = i [phi,psi]. With the tree routine below its f equals -field/2.
    # These substitutions retain the physical trace and all three coefficients.
    first = defaultdict(int)
    for (i,j,k),a,b in product(perms,axes,axes):
        add_to(first, mult(scalar_fermion[a][i],scalar_fermion[b][j],
                           fermion_triple[a][b][k]), eps(i,j,k))
    ff = pair(field,field)
    result = defaultdict(int)
    add_to(result, mult(*([ff] * n), clean(first)), Fraction(1, 4**n))
    if n:
        middle = defaultdict(int)
        last = defaultdict(int)
        for x,y in product(perms,repeat=2):
            i,j,k=x
            l,m,h=y
            sign=eps(*x)*eps(*y)
            for a in axes:
                add_to(middle, mult(field_fermion[i],scalar_fermion[a][l],
                                    fermion_pair[j][m],fermion_triple[a][k][h]), sign)
            for z in perms:
                add_to(last, mult(*(fermion_triple[a][b][c] for a,b,c in zip(x,y,z))),
                       sign*eps(*z))
        add_to(result, mult(*([ff]*(n-1)),clean(middle)), Fraction(n,4**n))
        add_to(result, mult(*([ff]*(n-1)),clean(last)), Fraction(2*n*n+n,432*4**(n-1)))
    return clean(result)


@cache
def hair(i):
    result = defaultdict(int)
    add_to(result, pair(phi[i],field), Fraction(-1,2))
    for j,k in product(axes,repeat=2):
        add_to(result, fermion_pair[j][k], Fraction(-eps(i,j,k),4))
    return clean(result)


def bracket(x,y,mu):
    result=defaultdict(int)
    for (a,b),(z,n) in frame.BRACKET.items():
        if mu==z:
            add_to(result,multiply(x[a],y[b]),n)
    return clean(result)


@cache
def replacement(op,axis,is_odd,pos):
    species,colour=divmod(pos,3)
    result=defaultdict(int)
    if op=='Q':
        if is_odd:
            return bracket(phi[(species+1)%3],phi[(species+2)%3],colour)
        if species==3:
            for i in axes:
                add_to(result,bracket(phi[i],psi[i],colour))
    elif op=='D':
        if is_odd and axis==species:
            return field[colour]
        if not is_odd and species!=3:
            for k in axes:
                add_to(result,psi[k][colour],eps(axis,species,k))
    elif op=='S':
        if is_odd:
            for k in axes:
                add_to(result,phi[k][colour],eps(axis,species,k))
        if not is_odd and species==3:
            add_to(result,psi[axis][colour],-2)
    else:
        raise ValueError(op)
    return clean(result)


def act(op,axis,poly):
    result=defaultdict(int)
    for (powers,mask),coefficient in poly.items():
        for pos in range(12):
            bit=1<<(4*pos)
            count=(powers//bit)&15
            if count:
                for (delta,newmask),n in replacement(op,axis,False,pos).items():
                    result[powers-bit+delta,mask|newmask] += (
                        coefficient*count*n*frame.wedge(newmask,mask))
        for pos in range(9):
            bit=1<<pos
            if mask&bit:
                sign=(-1)**((mask&(bit-1)).bit_count())
                for (delta,newmask),n in replacement(op,axis,True,pos).items():
                    assert newmask==0
                    result[powers+delta,mask^bit]+=sign*n*coefficient
    return clean(result)


def dressed(n):
    out=defaultdict(int)
    for i,j,k in perms:
        add_to(out,multiply(hair(k),act('D',i,act('D',j,core(n)))),eps(i,j,k))
    return clean(out)


def slice_vector(poly, auxiliary=(0,0,0)):
    out=defaultdict(int)
    for (powers,mask),coefficient in poly.items():
        if any((powers>>(4*i))&15 for i in (1,2,3,5,6,7)):
            continue
        f=tuple((powers>>(4*(9+i)))&15 for i in axes)
        out[mask,auxiliary,f] += coefficient
    return clean(out)


def integer_vector(row):
    den=lcm(*(getattr(c,'denominator',1) for c in row.values()))
    return {k:int(c*den) for k,c in row.items()}


def covariant(op,poly):
    out=defaultdict(int)
    for i in axes:
        add_to(out,slice_vector(act(op,i,poly),frame.unit(i)))
    return clean(out)


def boundaries(w,y):
    kind=frame.source_kind(w)
    d=sum(w)
    result=frame.Rank(0)
    for fdegree in range(d//3+1):
        k=d-y-2*fdegree
        phi_count=y-k-fdegree
        if not (0<=k<=9 and phi_count>=0):
            continue
        coeff=frame.coefficients(phi_count,k,kind,0)
        for f in frame.powers(fdegree):
            for v in coeff.get(-frame.weight(f),()):
                result.insert(frame.image(v,f))
    return result


def residue(basis, vec):
    v=integer_vector(vec)
    # Full reduction uses all pivots even when the first occupied column is free.
    for col,pivot in sorted(basis.rows.items()):
        c=v.get(col,0)
        if c:
            a=pivot[col]
            g=gcd(a,c)
            updated=defaultdict(int,{k:x*(a//g) for k,x in v.items()})
            add_to(updated,pivot,-(c//g))
            v=clean(updated)
            if v:
                div=gcd(*v.values())
                v={k:x//div for k,x in v.items()}
    return v


def joined(w, items, y=11):
    basis=boundaries(w,y-1)
    size=len(basis)
    quotient=frame.Rank(0)
    result={}
    for name,row in items.items():
        reduced=residue(basis,row)
        inc=quotient.insert(reduced)
        result[name]=dict(nonzero=bool(reduced),adds_dimension=inc,support=len(reduced))
    return dict(diagonal_weight_zero_incoming=size,additional_rank=len(quotient),vectors=result)


def main():
    for n in range(3):
        assert not act('Q',0,core(n)), n
    for i in axes:
        assert not act('Q',0,hair(i)) and all(not act('D',j,hair(i)) for j in axes)
    source={'O2':core(2),'epsilon_w_DD_O1':dressed(1)}
    assert all(not act('Q',0,v) for v in source.values())
    result={
        'singlet48':joined((8,8,8),{k:slice_vector(v) for k,v in source.items()}),
        'lower46':joined((8,8,7),{k:covariant('S',v) for k,v in source.items()}),
        'raise50':joined((9,8,8),{k:covariant('D',v) for k,v in source.items()}),
        'terms':{k:len(v) for k,v in source.items()}
    }
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    if not __debug__:
        raise ValueError('Assertions must be enabled')
    main()
