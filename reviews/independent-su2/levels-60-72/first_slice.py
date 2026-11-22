"""Exact first-image tests after setting the undotted scalars to zero and f=H.

One dotted derivative in each direction is projected to the spin singlet.
All possible tree boundaries on this specialization are included before
testing the projected ordered-pair kernel. Both maps use divided Taylor jets.
"""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, product
import json
import descendants as d


axes=range(3)
flavours=range(8)
colours=range(3)
letters=tuple(product(flavours,colours))


def parity(mask):
    return (mask.bit_count()+1)%2


def sign(order):
    return (-1)**sum(x>y for i,x in enumerate(order) for y in order[i+1:])


def components(s):
    return [i for i in axes if s & (1<<i)]


def perm_sign(current, target, odd):
    pos={v:i for i,v in enumerate(target)}
    return sign([pos[t] for t in current if t in odd])


@cache
def pair_rule(s,t):
    # Coefficient of theta^s eta^t in (theta-eta)^3 dC(theta)dC(eta).
    # The two dots are contracted with the alternating tensor.
    u=components(s)
    v=[3+i for i in components(t)]
    odd=set(range(6))|({6} if parity(s) else set())|({7} if parity(t) else set())
    left=(-1)**(s.bit_count()+t.bit_count())*perm_sign(u+[6]+v+[7],u+v+[6,7],odd)
    terms=[]
    for first in flavours:
        second=7^first
        if first & ~s or second & ~t:
            continue
        a=s^first
        b=t^second
        prefix=[i if first&(1<<i) else i+3 for i in axes]
        cur=prefix+components(a)+[6]+[3+i for i in components(b)]+[7]
        odd=set(range(6))|({6} if parity(a) else set())|({7} if parity(b) else set())
        n=left*(-1)**second.bit_count()*perm_sign(cur,u+v+[6,7],odd)
        terms.extend(((n,a,b,True),(-n,a,b,False)))
    return tuple(terms)


@cache
def colour_kernel(first,second):
    result=defaultdict(Fraction)
    for u,v,c in ((0,2,1),(1,1,Fraction(1,2)),(2,0,1)):
        for a,b in product(colours,repeat=2):
            x=d.frame.BRACKET.get((u,a))
            y=d.frame.BRACKET.get((v,b))
            if x and y and (x[0],y[0])==(first,second):
                result[a,b]+=x[1]*y[1]*c
    return d.clean(result)


def order_dots(a,b):
    spin_sign=-(-1)**(parity(a[0])*parity(b[0]))
    if a==b and spin_sign==-1:
        return 0,a,b
    return (1,a,b) if a<=b else (spin_sign,b,a)


@cache
def tree_factor(a,field_colour):
    # Retain the ordered product (dotted letter)(undotted letter) in Q C=C^2.
    # Its partner in the opposite order supplies the colour commutator.
    s,mu=a
    terms=[]
    for t in flavours:
        if t & ~s:
            continue
        u=s^t
        if u.bit_count()<2:
            continue
        c=(-1)**(s.bit_count()+parity(t)*u.bit_count())*sign(components(t)+components(u))
        # In the frame C_(jk) = eps_(ijk) psi_i, C_(123) = -f.
        i=(7^u).bit_length()-1 if u!=7 else None
        bridge=(-1)**i if i is not None else -1
        for (col,other),(dest,coef) in d.frame.BRACKET.items():
            if dest!=mu or (i is None and other!=field_colour):
                continue
            ext=1<<(3*i+other) if i is not None else 0
            terms.append((c*bridge*coef,(t,col),ext))
    return terms


def tree(key,field_colour=1):
    a,b,mask=key
    result=defaultdict(int)
    for n,new,e in tree_factor(a,field_colour):
        ss,x,y=order_dots(new,b)
        cross=(-1)**(int(bool(e))*parity(b[0]))
        result[x,y,mask|e]+=n*ss*cross*d.frame.wedge(e,mask)
    for n,new,e in tree_factor(b,field_colour):
        ss,x,y=order_dots(a,new)
        result[x,y,mask|e]+=(-1)**parity(a[0])*n*ss*d.frame.wedge(e,mask)
    return d.clean(result)


def basis(weight,y,field_colour=1):
    for j,a in enumerate(letters):
        for b in letters[j:]:
            if not order_dots(a,b)[0]:
                continue
            counts=[y-2+((a[0]>>i)&1)+((b[0]>>i)&1)-weight[i] for i in axes]
            if min(counts)<0 or max(counts)>3 or sum(counts)>y-2:
                continue
            fields=y-2-sum(counts)
            if fields and field_colour is None:
                continue
            for chosen in product(*(combinations(colours,n) for n in counts)):
                colour=d.frame.W[a[1]]+d.frame.W[b[1]]+sum(d.frame.W[c] for row in chosen for c in row)
                if colour+(fields*d.frame.W[field_colour] if fields else 0):
                    continue
                mask=sum(1<<(3*i+c) for i,row in enumerate(chosen) for c in row)
                yield a,b,mask


def as_letters(powers,mask):
    out=[]
    for pos in range(12):
        n=(powers>>(4*pos))&15
        typ,col=divmod(pos,3)
        out.extend([((1<<typ if typ<3 else 7),col)]*n)
    for i in range(9):
        if mask&(1<<i):
            typ,col=divmod(i,3)
            out.append((7^(1<<typ),col))
    return out


def quantum(poly,field_colour=1):
    result=defaultdict(Fraction)
    for (powers,mask),coefficient in poly.items():
        factors=as_letters(powers,mask)
        must_remove={j for j,(s,c) in enumerate(factors)
                     if s.bit_count()==1 or (s==7 and c!=field_colour)}
        if len(must_remove)>2:
            continue
        for i,j in combinations(range(len(factors)),2):
            if not must_remove <= {i,j}:
                continue
            (s,mu),(t,nu)=factors[i],factors[j]
            c=coefficient
            residual=mask
            for pos in (i,j):
                typ,col=factors[pos]
                if typ==7:
                    c=-c
                elif typ.bit_count()==2:
                    h=(7^typ).bit_length()-1
                    c*=(-1)**h
                    bit=1<<(3*h+col)
                    c*=(-1)**((residual & (bit-1)).bit_count())
                    residual^=bit
            for n,a,b,plus_first in pair_rule(s,t):
                for (ca,cb),z in colour_kernel(mu,nu).items():
                    x,y=(a,ca),(b,cb)
                    if plus_first:
                        k,one,two=order_dots(x,y)
                    else:
                        k,one,two=order_dots(y,x)
                        k*=(-1)**(parity(a)*parity(b))
                    result[one,two,residual]+=k*n*c*z
    return d.clean(result)


def evaluate(n,field_colour=1,hairy=False):
    y=7+2*n
    weight=(3+2*n,)*3
    source={'O'+str(n):d.core(n)}
    if n:
        source['epsilon_w_DD_O'+str(n-1)]=d.dressed(n-1)
    if n==2:
        primary=defaultdict(int)
        d.add_to(primary,d.core(2),102)
        d.add_to(primary,source['epsilon_w_DD_O1'],-1)
        source['102_O2_minus_epsilon_w_DD_O1']=d.clean(primary)
    if hairy:
        source={'w_'+name:d.multiply(d.hair(0),pol) for name,pol in source.items()}
        y+=2
        weight=tuple(x+1+(i==0) for i,x in enumerate(weight))
    target=set(basis(weight,y,field_colour))
    boundary=d.frame.Rank(0)
    before=list(basis(weight,y-1,field_colour))
    for x in before:
        v=tree(x,field_colour)
        assert set(v)<=target
        boundary.insert(v)
    output=d.frame.Rank(0)
    images={}
    for name,term in source.items():
        v=quantum(term,field_colour)
        assert set(v)<=target,(name,set(v)-target)
        after=defaultdict(int)
        for x,c in v.items():
            d.add_to(after,tree(x,field_colour),c)
        assert not d.clean(after),name
        r=d.residue(boundary,v)
        inc=output.insert(r)
        images[name]=dict(coordinates=len(v),residual_support=len(r),adds_dimension=inc)
    return dict(weight=weight,predecessors=len(before),target=len(target),
                rational_boundary_rank=len(boundary),rational_detected_rank=len(output),images=images)


if __name__=='__main__':
    if not __debug__:
        raise ValueError('Assertions must be enabled')
    print(json.dumps({n:evaluate(n) for n in range(3)},indent=2))
