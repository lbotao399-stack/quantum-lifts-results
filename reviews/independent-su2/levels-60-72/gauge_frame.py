"""Complete flavour and colour singlet calculation on two dense frames."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations_with_replacement
from math import gcd, lcm
from pathlib import Path
import json
import sys
from time import monotonic

import first_frame as frame
from first_frame import d,native


@cache
def changes(generator, slot, dual):
    bracket=d.frame.BRACKET
    if dual:
        return tuple((new,-c) for (gen,old),(new,c) in bracket.items()
                     if gen==generator and old==slot)
    return tuple((old,c) for (gen,old),(new,c) in bracket.items()
                 if gen==generator and new==slot)


@cache
def change_jet(generator,jet):
    s,mu=jet
    out=[((s,v),c) for v,c in changes(generator,mu,False)]
    if s.bit_count()==1:
        i=s.bit_length()-1
        out.extend(((1<<v,mu),c) for v,c in changes(generator,i,True))
    elif s.bit_count()==2:
        i=(7^s).bit_length()-1
        out.extend(((7^(1<<v),mu),c*(-1)**(i+v)) for v,c in changes(generator,i,False))
    return out


@cache
def colour_column(generator,a,b,mask):
    out=defaultdict(int)
    for position,jet in enumerate((a,b)):
        for new,c in change_jet(generator,jet):
            sign,x,y=native.order_dots(a,new) if position else native.order_dots(new,b)
            out[x,y,mask]+=sign*c
    for pos in range(9):
        bit=1<<pos
        if not mask&bit:continue
        i,mu=divmod(pos,3)
        remaining=mask^bit
        contributions=[(3*i+v,c) for v,c in changes(generator,mu,False)]
        contributions += [(3*v+mu,c) for v,c in changes(generator,i,False)]
        for index,c in contributions:
            new=1<<index
            phase=(-1)**((remaining&(bit-1)).bit_count())
            out[a,b,remaining|new]+=phase*d.frame.wedge(new,remaining)*c
    return d.clean(out)


def colour(row, gen):
    result=defaultdict(int)
    for key,c in row.items():
        d.add_to(result,colour_column(gen,*key),c)
    return d.clean(result)


@cache
def casimir_column(a,b,mask):
    return colour(colour_column(0,a,b,mask),2)


def primitive(row):
    if not row:return {}
    den=lcm(*(getattr(v,'denominator',1) for v in row.values()))
    integral={k:int(v*den) for k,v in row.items() if v}
    content=gcd(*integral.values())
    return {k:v//content for k,v in integral.items()}


def invariant(k,m,i,j,nf,prime):
    rows=frame.cochain_rows(k,m,i,j,prime)
    origin=rows.get(0,())
    if not origin:return ()
    highest=max(rows)
    multiplicities=[len(rows.get(v,()))-len(rows.get(v+1,())) for v in range(highest+1)]
    assert min(multiplicities)>=0
    allowed=sum(n for ell,n in enumerate(multiplicities) if ell<=nf and (nf-ell)%2==0)
    if not allowed:return ()
    for ell in reversed(range(highest+1)):
        if not multiplicities[ell] or (ell<=nf and (nf-ell)%2==0):
            continue
        rank=d.frame.Rank(prime)
        for row in origin:
            image=defaultdict(int)
            for key,c in row.items():
                d.add_to(image,casimir_column(*key),c)
                image[key]-=ell*(ell+1)*c
            rank.insert(image)
        origin=tuple(rank.rows.values())
    assert len(origin)==allowed,(k,m,i,j,nf,multiplicities,len(origin),allowed)
    return origin


def profiles(w,y):
    for i,j in combinations_with_replacement(range(4),2):
        for nf in range(y-1):
            k=3*w-(y-2)-2*nf-i-j
            m=y-2-k-nf
            if 0<=k<=9 and m>=0:
                yield k,m,i,j,nf


@cache
def tree_at_fixed_field(a,b,mask):
    result=defaultdict(int)
    for (x,y,ext,powers),c in frame.tree_column(a,b,mask,(0,0,0)).items():
        if powers[0]==powers[2]==0:
            result[x,y,ext]+=c
    return d.clean(result)


def differential_twice(row,nf):
    # Clear the two transverse denominators before accumulating exact entries.
    result=defaultdict(int)
    for key,c in row.items():
        d.add_to(result,tree_at_fixed_field(*key),2*c)
    # [E,H]=-2E and [F,H]=2F. Homogeneity determines the radial derivative.
    for mu,derivative,scale in ((0,colour(row,0),1),(1,row,2*nf),
                                 (2,colour(row,2),-1)):
        if not scale:continue
        for (a,b,mask),coef in derivative.items():
            phase=(-1)**(native.parity(a[0])+native.parity(b[0]))
            for pos,n in d.frame.Q_F[mu].items():
                bit=1<<pos
                result[a,b,mask|bit]+=phase*scale*coef*n*d.frame.wedge(bit,mask)
    return d.clean(result)


def full_tree(row,w,y):
    groups=defaultdict(dict)
    for key,c in row.items():
        a,b,mask=key
        count=3*w-(y-2)-a[0].bit_count()-b[0].bit_count()-mask.bit_count()
        assert count>=0 and count%2==0
        groups[count//2][key]=c
    result=defaultdict(int)
    for nf,terms in groups.items():
        d.add_to(result,differential_twice(terms,nf))
    return d.clean(result)


def scalar_field_slice(row):
    result=defaultdict(Fraction)
    for (a,b,mask,powers),c in row.items():
        if powers[0]==powers[2]==0:
            result[a,b,mask]+=c
    return d.clean(result)


def calculate(n,prime,check_all):
    start=monotonic()
    y,w=7+2*n,3+2*n
    boundary=d.frame.Rank(prime)
    tested=0
    test=lambda row: all(int(v.numerator)%prime==0 for v in row.values()) if prime else not row
    for k,m,i,j,nf in profiles(w,y-1):
        basis=invariant(k,m,i,j,nf,prime)
        for idx,vec in enumerate(basis):
            image=differential_twice(vec,nf)
            if not idx or check_all:
                assert test(full_tree(image,w,y)),(k,m,i,j,nf,idx)
            boundary.insert(d.integer_vector(image))
        tested+=len(basis)
        if basis:
            print(i,j,nf,len(basis),tested,len(boundary),round(monotonic()-start),
                  file=sys.stderr,flush=True)
    size=len(boundary)
    expected={1:584,2:1035,3:1195,4:1210}[n]
    assert tested==expected
    first=[]
    for src in (d.core(n),d.dressed(n-1)):
        assert not d.act('Q',0,src)
        image=scalar_field_slice(frame.quantum(src))
        assert not full_tree(image,w,y)
        first.append(dict(support=len(image),adds_rank=boundary.insert(d.integer_vector(image))))
    return dict(level=24+12*n,prime=prime,source_cochains=tested,
                incoming_tree_rank=size,joined_rank=len(boundary),images=first,
                first_rank=len(boundary)-size)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('n',type=int,choices=range(1,5))
    p.add_argument('--prime',type=int,default=1000003)
    p.add_argument('--check-all',action='store_true')
    q=p.parse_args()
    if not __debug__:p.error('assertions required')
    print(json.dumps(calculate(q.n,q.prime,q.check_all),indent=2))
