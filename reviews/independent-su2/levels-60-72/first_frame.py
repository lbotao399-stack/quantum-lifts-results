"""First supercharge on the complete dotted-singlet target, in a scalar frame."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path
import sys
from time import monotonic

import descendants as d
import first_slice as native


@cache
def masks(size, count):
    if not 0 <= count <= size:
        return ()
    return tuple(sum(1 << i for i in choice) for choice in combinations(range(size), count))


@cache
def determinant(kinds):
    terms = defaultdict(int)
    # Each of the five kinds has three odd, covariant flavour components.
    for perm in permutations(range(3)):
        sign, current = native.sign(perm), 0
        for kind, flavour in zip(kinds, perm):
            bit = 1 << (3*kind + flavour)
            sign *= d.frame.wedge(current, bit)
            current |= bit
        terms[current] += sign
    return d.clean(terms)


@cache
def flavour_rows(k, left, right, determinants, prime):
    if not 0 <= k <= 9 or min(left, right) < 0:
        return ()
    if not determinants:
        return tuple({i | (a << 9) | (b << 12): 1}
                     for i,a,b in product(masks(9,k), masks(3,left), masks(3,right)))
    rank = d.frame.Rank(prime)
    for kinds in combinations_with_replacement(range(5),3):
        a,b = kinds.count(3),kinds.count(4)
        ordinary = 3-a-b
        if a>left or b>right or ordinary>k:
            continue
        for row in flavour_rows(k-ordinary, left-a, right-b, determinants-1, prime):
            varied = defaultdict(int)
            for u,x in determinant(kinds).items():
                for v,y in row.items():
                    if not u & v:
                        varied[u|v] += d.frame.wedge(u,v)*x*y
            rank.insert(varied)
    return tuple(rank.rows.values())


def dotted(index, selector):
    if index in (0,3):
        assert selector == 0
        return (0 if index==0 else 7),1
    flavour_bit = selector if index==2 else 7^selector
    assert flavour_bit.bit_count()==1
    flavour = flavour_bit.bit_length()-1
    return (7^flavour_bit if index==2 else flavour_bit),(-1)**flavour


@cache
def coordinate(word, i, j, colour_a, colour_b):
    a,sa=dotted(i,(word>>9)&7)
    b,sb=dotted(j,(word>>12)&7)
    sign, first, second=native.order_dots((a,colour_a),(b,colour_b))
    return (first,second,ordinary_native(word&511)),sign*sa*sb


@cache
def ordinary_native(encoded):
    # Auxiliary basis is ordered by colour, then flavour. Physical exterior
    # variables are ordered by flavour, then colour.
    bits=[3*flavour+colour for colour in range(3) for flavour in range(3)
          if encoded & (1<<(3*colour+flavour))]
    return sum(1<<v for v in bits),native.sign(bits)


@cache
def joint_weight(key):
    a,b,mask=key
    value=d.frame.maskweight(mask)
    for s,col in (a,b):
        value += d.frame.W[col]
        if s.bit_count()==1:
            value -= d.frame.W[s.bit_length()-1]
        elif s.bit_count()==2:
            value += d.frame.W[(7^s).bit_length()-1]
    return value


@cache
def cochain_rows(k, m, i, j, prime):
    l,r=[(0,2,1,0)[q] for q in (i,j)]
    surplus=k+l+r-m
    assert surplus % 3 == 0
    base=flavour_rows(k,l,r,max(0,surplus//3),prime)
    rank=defaultdict(lambda:d.frame.Rank(prime))
    for ca,cb in product(range(3),repeat=2):
        if i==j and ca>cb:
            continue
        for row in base:
            values=defaultdict(int)
            for word,n in row.items():
                (a,b,(mask,c)),sg=coordinate(word,i,j,ca,cb)
                if sg:
                    values[a,b,mask]+=n*sg*c
            values=d.clean(values)
            if values:
                weights={joint_weight(x) for x in values}
                assert len(weights)==1
                rank[weights.pop()].insert(values)
    return {w:tuple(v.rows.values()) for w,v in rank.items()}


def basis(w,y,prime,charge):
    for i,j in combinations_with_replacement(range(4),2):
        for nf in range(y-1):
            k=3*w-(y-2)-2*nf-i-j
            m=y-2-k-nf
            if not 0<=k<=9 or m<0:
                continue
            rows=cochain_rows(k,m,i,j,prime)
            for powers in d.frame.powers(nf):
                for row in rows.get(charge-d.frame.weight(powers),()):
                    yield row,powers


@cache
def dotted_split(letter):
    s,dest=letter
    result=[]
    for t in range(8):
        if t & ~s or t==s:
            continue
        u=s^t
        sign=(-1)**(s.bit_count()+native.parity(t)*u.bit_count())
        sign*=native.sign(native.components(t)+native.components(u))
        for (col,other),(out,coef) in d.frame.BRACKET.items():
            if out!=dest:
                continue
            if u.bit_count()==1:
                if other==u.bit_length()-1: result.append((sign*coef,(t,col),0,-1))
            elif u.bit_count()==2:
                flavour=(7^u).bit_length()-1
                result.append((sign*coef*(-1)**flavour,(t,col),1<<(3*flavour+other),-1))
            else:
                result.append((-sign*coef,(t,col),0,other))
    return tuple(result)


@cache
def tree_column(a,b,mask,fields):
    values=defaultdict(int)
    for pos,item in enumerate((a,b)):
        for c,changed,ext,col in dotted_split(item):
            orientation,x,y=native.order_dots(a,changed) if pos else native.order_dots(changed,b)
            phase=(-1)**(native.parity(a[0]) if pos else bool(ext)*native.parity(b[0]))
            powers=fields if col==-1 else tuple(t+(j==col) for j,t in enumerate(fields))
            values[x,y,mask|ext,powers]+=c*phase*orientation*d.frame.wedge(ext,mask)
    phase=(-1)**(native.parity(a[0])+native.parity(b[0]))
    for (changed,powers),c in d.frame.on_monomial(mask,fields).items():
        values[a,b,changed,powers]+=phase*c
    return d.clean(values)


def tree(row):
    values=defaultdict(int)
    for key,c in row.items():
        d.add_to(values,tree_column(*key),c)
    return d.clean(values)


def quantum(poly):
    values=defaultdict(Fraction)
    for (encoded,mask),c in poly.items():
        letters=native.as_letters(encoded,mask)
        required={i for i,(s,col) in enumerate(letters)
                  if s.bit_count()==1 and s.bit_length()-1!=col}
        if len(required)>2:
            continue
        fields=tuple((encoded>>(4*(9+i)))&15 for i in range(3))
        for i,j in combinations(range(len(letters)),2):
            if not required <= {i,j}:
                continue
            (s,mu),(t,nu)=letters[i],letters[j]
            residual,powers,value=mask,fields,c
            for v,col in (letters[i],letters[j]):
                if v==7:
                    value=-value
                    powers=tuple(n-(k==col) for k,n in enumerate(powers))
                elif v.bit_count()==2:
                    h=(7^v).bit_length()-1
                    bit=1<<(3*h+col)
                    value*=(-1)**(h+((residual & (bit-1)).bit_count()))
                    residual^=bit
            for n,u,v,plus in native.pair_rule(s,t):
                for (cu,cv),z in native.colour_kernel(mu,nu).items():
                    one,two=(u,cu),(v,cv)
                    sign,x,y=native.order_dots(one,two) if plus else native.order_dots(two,one)
                    if not plus: sign*=(-1)**(native.parity(u)*native.parity(v))
                    values[x,y,residual,powers]+=n*value*sign*z
    return d.clean(values)


def calculate(n, prime, dimension_only):
    w=3+2*n
    y=7+2*n
    start=monotonic()
    dims=[sum(1 for _ in basis(w,y-1,prime,z)) for z in (0,1)]
    print('predecessor dimensions',dims,'time',round(monotonic()-start),file=sys.stderr,flush=True)
    if dimension_only:
        return dict(level=24+12*n,prime=prime,predecessor_weight_dimensions=dims)
    boundary=d.frame.Rank(prime)
    for idx,(row,powers) in enumerate(basis(w,y-1,prime,0)):
        vector={(*key,powers):c for key,c in row.items()}
        boundary.insert(tree(vector))
        if (idx+1)%2500==0:
            print(idx+1,'boundary rank',len(boundary),round(monotonic()-start),file=sys.stderr,flush=True)
    assert idx+1==dims[0]
    before=len(boundary)
    supports=[]
    for src in (d.core(n),d.dressed(n-1)):
        assert not d.act('Q',0,src)
        image=quantum(src)
        assert not tree(image)
        supports.append(len(image))
        boundary.insert(d.integer_vector(image))
    return dict(level=24+12*n,prime=prime,predecessor_weight_dimensions=dims,
                boundary_rank=before,first_image_dimensions=len(boundary)-before,
                first_image_supports=supports)


if __name__=='__main__':
    a=ArgumentParser(description=__doc__)
    a.add_argument('n',type=int,choices=range(1,5))
    a.add_argument('--prime',type=int,default=1000003)
    a.add_argument('--dimension-only',action='store_true')
    q=a.parse_args()
    if not __debug__:a.error('assertions required')
    print(json.dumps(calculate(q.n,q.prime,q.dimension_only),indent=2))
