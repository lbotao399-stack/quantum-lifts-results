"""Classify the first syzygies of one F4 commutator with one scalar of each flavour."""
from collections import Counter
from itertools import combinations, product
import json

from colour import Q, add, compact, items
from second_test import Second

def weights_sum(*weights):
    return tuple(map(sum,zip(*weights)))

def orbits(colour):
    multiplicities = Counter(weights_sum(*ws) for ws in product(colour.weights,repeat=3))
    coroots = []
    for e,f in zip(colour.data['raising'],colour.data['lowering']):
        root = colour.weights[e]
        h = colour.br[e][f]
        denom = sum(root[i]*x for i,x in h.items())
        coroots.append((tuple(Q(2*h.get(i,0),denom) for i in range(4)),root))
    unseen = set(multiplicities)
    while unseen:
        start = max(unseen)
        orbit,pending = {start},[start]
        while pending:
            w = pending.pop()
            for h,root in coroots:
                n = sum(x*y for x,y in zip(h,w))
                assert Q(n).denominator==1
                image = tuple(x-int(n)*y for x,y in zip(w,root))
                if image not in orbit:
                    assert multiplicities[image]==multiplicities[start]
                    orbit.add(image)
                    pending.append(image)
        dominant = [w for w in orbit if all(sum(x*y for x,y in zip(h,w))>=0 for h,_ in coroots)]
        assert len(dominant)==1 and orbit<=unseen
        unseen-=orbit
        yield dominant[0],len(orbit)

def known_cycles(colour,weight,cols):
    idx = {v:i for i,v in enumerate(cols)}
    koszul,products = [],[]
    for a,b in combinations(range(52),2):
        if weights_sum(colour.weights[a],colour.weights[b])!=weight:
            continue
        row = {}
        for r,s,factor in ((a,b,1),(b,a,-1)):
            for mask,x in colour.delta[s].items():
                d,e = items(mask)
                add(row,idx[r,d,e],factor*x)
                add(row,idx[r,e,d],-factor*x)
        koszul.append(compact(row))
    for r,w in enumerate(colour.weights):
        if w==weight:
            products.append({idx[a,d,r]:x for (a,d),x in colour.inverse.items()})
            products.append({idx[a,r,d]:x for (a,d),x in colour.inverse.items()})
    return koszul,products

def rank(rows):
    pivots = {}
    for vec in rows:
        row = vec.copy()
        while row:
            i = min(row)
            x = row[i]
            if i not in pivots:
                pivots[i] = {k:Q(v)/x for k,v in row.items()}
                break
            for k,v in pivots[i].items():
                add(row,k,-x*v)
                if not row[k]:
                    del row[k]
    return len(pivots)

def classify():
    engine = Second()
    colour = engine.c
    blocks = []
    for weight,size in orbits(colour):
        cols,rows,kernel = engine.classical(weight)
        boundaries,products = known_cycles(colour,weight,cols)
        known = boundaries+products
        # The computed kernel was checked column by column by classical().
        # Equality of these two spans checks closure as well as completeness.
        k,r,j = rank(kernel),rank(known),rank(kernel+known)
        assert k == len(kernel) and r == len(known) and j == k
        blocks.append(dict(weight=weight,orbit=size,columns=len(cols),output_rows=rows,
                           kernel=k,koszul=len(boundaries),bilinear_products=len(products),unexplained=k-r))
    fields = ('columns','kernel','koszul','bilinear_products','unexplained')
    total = {k:sum(b['orbit']*b[k] for b in blocks) for k in fields}
    assert total['columns']==52**3 and total['koszul']==52*51//2 and total['bilinear_products']==2*52
    return dict(weyl_orbits=sorted(blocks,key=lambda b:b['weight']),total=total)
if __name__=='__main__':
    if not __debug__:
        raise SystemExit('Run with assertions enabled')
    print(json.dumps(classify(),indent=2))
