"""Bound homology in one Taylor and colour representation by two augmented ranks."""
from collections import Counter
from functools import cache
from math import comb
from determinant import basis, boundary, packet

def dominant(w):
    a,b=w
    while min(a,b)<0:
        if a<0:a,b=-a,b+a
        else:a,b=a+3*b,-b
    return a,b

def weyl():
    signs={(1,1):1}
    new=[(1,1)]
    while new:
        a,b=new.pop()
        for w in ((-a,a+b),(a+3*b,-b)):
            if w not in signs:
                signs[w]=-signs[a,b]
                new.append(w)
    return signs

def determinant_dimensions(ce, other, length, weight):
    @cache
    def amount(sizes,x,y):
        if not sizes:return int(x==y==0)
        if len(sizes)==1:return len(ce.exterior(sizes[0]).get((x,y),()))
        return sum(len(group)*amount(sizes[1:],x-a,y-b) for (a,b),group in ce.exterior(sizes[0]).items())
    def raw(pair,w):
        return sum(amount(tuple(sorted((n for m,n in pat),key=lambda k:comb(14,k))),*w)
                   for pat in pair.patterns(length))
    x,y=weight
    dims={}
    for (a,b),sgn in weyl().items():
        w=dominant((x+1-a,y+1-b))
        dims[w]=raw(ce,w)-raw(other,w)
    return dims[weight],sum(sgn*dims[dominant((x+1-a,y+1-b))] for (a,b),sgn in weyl().items())

def lower(words, projection, algebra, root):
    colour=algebra.weights.index([-n for n in root])
    for word in words:
        col=Counter()
        for k,(d,c,a) in enumerate(word):
            for g,n in algebra.brackets[colour][c]:
                changed=(d,g,a)
                rest=word[:k]+word[k+1:]
                if changed in rest:continue
                groups,state,sign=packet(word[:k]+(changed,)+word[k+1:])
                offset,mapping=projection[groups]
                for j,v in mapping[state].items():
                    col[offset+j]+=n*sign*v
        yield {i:n for i,n in col.items() if n}

def calculate(ce, other, weight, primes):
    from panels import rank
    l=ce.length
    dims=[determinant_dimensions(ce,other,k,weight) for k in (l-1,l,l+1)]
    spaces=[basis(ce,k,weight) for k in (l-1,l,l+1)]
    assert [len(words) for words,_ in spaces]==[d[0] for d in dims]
    maps=[]
    for i in (0,1):
        dst=spaces[i][1]
        cols=list(boundary(spaces[i+1][0],dst,ce.brackets))
        maps.append(cols)
    for col in maps[1]:
        square=Counter()
        for j,m in col.items():
            for i,n in maps[0][j].items(): square[i]+=m*n
        assert not any(square.values())
    out={}
    for i,(words,projection) in enumerate(spaces[:2]):
        k=l-1+i
        cols=maps[i]
        counts=[len(cols)]
        for root in ((2,-1),(-3,2)):
            x,y=weight
            source,_=basis(ce,k,(x+root[0],y+root[1]))
            counts.append(len(source))
            cols.extend(lower(source,projection,ce,root))
        qr=dims[i][0]-dims[i][1]
        r=[rank(cols,p) for p in primes]
        bound=max(0,max(r)-qr)
        assert bound<=min(dims[i][1],dims[i+1][1])
        out[k]=dict(column_counts=counts,relation_dimension=qr,augmented_ranks=r,rank_lower=bound)
    h=dims[1][1]-out[l-1]['rank_lower']-out[l]['rank_lower']
    assert h>=0
    return dict(weight=weight,primes=primes,chain_dimensions=dims,maps=out,rational_homology_upper=h)
