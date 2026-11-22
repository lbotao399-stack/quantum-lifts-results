"""Exact Cartan multiplication on the second local residual of the level-34 syzygy."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import gcd, lcm

from flint import fmpz_mat


def clean(row):
    return {w:c for w,c in row.items() if c}


def add(*pieces):
    result = Counter()
    for scale,row in pieces:
        for w,c in row.items():
            result[w] += scale*c
    return clean(result)


def odd(jet):
    return (1 + jet[2].bit_count()) % 2


def canonical(sequence):
    fermions = [j for j in sequence if odd(j)]
    if len(set(fermions)) != len(fermions):
        return (),0
    swaps = sum(x > y for i,x in enumerate(fermions) for y in fermions[i+1:])
    return tuple(sorted(sequence)),(-1)**swaps


def mul(left,right):
    result = Counter()
    for w,c in left.items():
        for v,d in right.items():
            key,sign = canonical(w+v)
            if sign:
                result[key] += c*d*sign
    return clean(result)


@lru_cache(None)
def trace(masks):
    polynomial = Counter()
    for colours in product(range(2),repeat=len(masks)):
        # X=diag(x,y,-x-y), with independent traceless coordinates for each jet.
        scale = int(not any(colours)) + int(all(colours)) + (-1)**len(masks)
        w,s = canonical(tuple((0,0,m,c) for m,c in zip(masks,colours,strict=True)))
        if s:
            polynomial[w] += s*scale
    return clean(polynomial)


def raw(upper,lower):
    # psi_1, psi_2, psi_3 have wedge phases +,-,+ in theta coordinates.
    masks = tuple(1 << i for i in upper) + (7 ^ (1 << lower),)
    return add(((-1)**lower,trace(masks)))


@lru_cache(None)
def v2(i,j):
    return add((1,raw((i,),j)), *[(Fraction(-1,3),raw((k,),k)) for k in range(3) if i == j])


@lru_cache(None)
def v3(i,j,k):
    return add((1,raw((i,j),k)), *[(Fraction(-1,4),raw((j,h),h)) for h in range(3) if i == k],
               *[(Fraction(-1,4),raw((i,h),h)) for h in range(3) if j == k])


def prefactors():
    # Highest free upper indices 0,0,0 and lower index 2. The paper orders v2 before v3.
    parts = [[] for _ in range(3)]
    for k in range(3):
        parts[0].append((199920,mul(v2(0,k),v3(k,0,2))))  # s37
        parts[k].append((-97104,mul(v2(0,k),v3(0,0,2))))  # s38
        parts[k].append((237048,mul(v2(0,2),v3(0,0,k))))   # s40
    return [add(*part) for part in parts]


def monomials():
    target = (2,2,2,1,1)
    letters = sorted((p,m,t,c) for p,m,t,c in product(range(3),range(3),range(8),range(2)) if p+m+t)
    weights = [(p,m,*((t>>i)&1 for i in range(3))) for p,m,t,c in letters]
    def visit(start,remaining,degree):
        if not degree:
            if not any(remaining):
                yield ()
            return
        if sum(remaining) < degree:
            return
        for i in range(start,len(letters)):
            rest = tuple(a-b for a,b in zip(remaining,weights[i]))
            if min(rest) >= 0:
                for suffix in visit(i+odd(letters[i]),rest,degree-1):
                    yield (letters[i],)+suffix
    return visit(0,target,5)


def profile(w):
    return tuple(sorted((p+m,t.bit_count()) for p,m,t,c in w))


def equations(w):
    out = [Counter() for _ in range(5)]
    # Simple Weyl reflections, followed by the two flavour roots and the dotted root.
    v,s = canonical(tuple((p,m,t,1-c) for p,m,t,c in w))
    out[0][v] += s
    out[0][w] -= 1
    rows = {():1}
    for p,m,t,c in w:
        if c:
            factor = {((p,m,t,0),):-1,((p,m,t,1),):-1}
        else:
            factor = {((p,m,t,0),):1}
        rows = mul(rows,factor)
    out[1].update(rows)
    out[1][w] -= 1
    for j,(p,m,t,c) in enumerate(w):
        for i in (0,1):
            if (t & (2<<i)) and not (t & (1<<i)):
                v,s = canonical(w[:j]+((p,m,t^(3<<i),c),)+w[j+1:])
                if s:
                    out[2+i][v] += s
        if m:
            v,s = canonical(w[:j]+((p+1,m-1,t,c),)+w[j+1:])
            if s:
                out[4][v] += s*(p+1)
    return {(i,v):c for i,row in enumerate(out) for v,c in row.items() if c}


def rotate(w,k):
    items,sign = [],1
    for p,m,t,c in w:
        indices = [(i+k)%3 for i in range(3) if t & (1<<i)]
        sign *= (-1)**sum(x>y for j,x in enumerate(indices) for y in indices[j+1:])
        items.append((p,m,sum(1<<i for i in indices),c))
    v,s = canonical(tuple(items))
    return v,s*sign


def matrix(columns):
    labels = sorted({w for col in columns for w in col})
    values = [int(col.get(w,0)) for w in labels for col in columns]
    return fmpz_mat(len(labels),len(columns),values)


def basis_for(columns):
    equations_matrix = matrix(columns)
    basis,count = equations_matrix.nullspace()
    for j in range(count):
        values = [int(basis[i,j]) for i in range(basis.nrows())]
        common = gcd(*values)
        yield [v//common for v in values]


def cartan_character():
    """Independent Molien and highest-weight extraction for five Cartan letters."""
    theta = (2,1,1)
    weights = Counter()
    for perm in permutations(range(3)):
        shifted = tuple(theta[i]+perm[i]-i for i in range(3))
        if min(shifted) < 0:
            continue
        sign = (-1)**sum(x > y for i,x in enumerate(perm) for y in perm[i+1:])
        weights[2,2,*shifted] += sign
        weights[3,1,*shifted] -= sign
    letters = sorted({(p,m,*t) for w in weights for p,m in product(range(w[0]+1),range(w[1]+1))
                      for t in product(range(2),repeat=3) if p+m+sum(t) and all(x<=y for x,y in zip(t,w[2:]))},
                     key=sum,reverse=True)

    @lru_cache(None)
    def multiplicities(i,n,w):
        if n == 0:
            return (int(not any(w)),)*3
        if i == len(letters) or sum(w) < n:
            return (0,0,0)
        letter = letters[i]
        fermion = (1+sum(letter[2:]))%2
        bound = min(n,2 if fermion else n,*(a//b for a,b in zip(w,letter) if b))
        result = [0,0,0]
        for k in range(bound+1):
            # S3 characters in the order identity, transposition, three-cycle.
            factor = ([(1,1,1),(2,0,-1),(1,-1,1)][k] if fermion else
                      (k+1,int(k%2 == 0),(1,-1,0)[k%3]))
            smaller = multiplicities(i+1,n-k,tuple(a-k*b for a,b in zip(w,letter)))
            for g in range(3):
                result[g] += factor[g]*smaller[g]
        return tuple(result)

    traces = [sum(sign*multiplicities(0,5,w)[g] for w,sign in weights.items()) for g in range(3)]
    numerator = traces[0] + 3*traces[1] + 2*traces[2]
    assert numerator % 6 == 0
    return traces,numerator//6


def calculate():
    factors = prefactors()
    by_profile = defaultdict(list)
    for w in monomials():
        by_profile[profile(w)].append(w)
    source_dims = image_ranks = 0
    noninjective = []
    for key,words in sorted(by_profile.items()):
        images = []
        for coefficients in basis_for([equations(w) for w in words]):
            source = {w:c for w,c in zip(words,coefficients,strict=True) if c}
            terms = []
            for k,multiplier in enumerate(factors):
                component = {}
                for w,c in source.items():
                    v,sign = rotate(w,k)
                    component[v] = sign*c
                terms.append((1,mul(multiplier,component)))
            image = add(*terms)
            expected = tuple(sorted(key + ((0,1),)*3 + ((0,2),)*2))
            assert all(profile(w) == expected for w in image)
            denominator = lcm(*(v.denominator for v in image.values()))
            images.append({w:int(c*denominator) for w,c in image.items()})
        if images:
            rank = matrix(images).rank()
            if rank < len(images):
                noninjective.append((key,len(images),rank))
            image_ranks += rank
            source_dims += len(images)
    traces,expected_dimension = cartan_character()
    assert source_dims == expected_dimension
    return dict(raw=sum(map(len,by_profile.values())),profiles=len(by_profile),weyl_traces=traces,
                prefactor_terms=[len(f) for f in factors],source_dimension=source_dims,
                image_rank=image_ranks,noninjective=noninjective)


if __name__ == '__main__':
    import json
    if not __debug__:
        raise SystemExit('enable assertions')
    result = calculate()
    assert (result['raw'],result['profiles'],result['prefactor_terms']) == (4678,20,[36,14,4])
    assert result['source_dimension'] == result['image_rank'] == 53 and not result['noninjective']
    print(json.dumps(result,indent=2))
