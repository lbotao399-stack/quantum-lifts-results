"""Certify a G2 exterior-algebra relation and its consequence for Taylor-cohomology lengths."""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import json
from math import comb


def structure_constants():
    triples = ((0,1,2,1),(0,3,4,1),(0,5,6,1),(1,3,5,1),
               (1,4,6,-1),(2,3,6,-1),(2,4,5,-1))
    form = defaultdict(int)
    for i,j,k,s in triples:
        for x,y,z in ((i,j,k),(j,k,i),(k,i,j)):
            form[x,y,z],form[y,x,z] = s,-s
    basis, coordinates = [], []
    for k in range(7):
        pairs = [(i,j) for i,j in combinations(range(7),2) if form[i,j,k]]
        assert len(pairs)==3
        for i,j in pairs[:2]:
            matrix = [[0]*7 for _ in range(7)]
            for x,y,sign in ((i,j,1),(*pairs[2],-1)):
                matrix[x][y] = sign*form[x,y,k]
                matrix[y][x] = -matrix[x][y]
            coordinates.append((i,j,matrix[i][j]))
            basis.append(matrix)
    assert len(basis)==14
    result = [[] for _ in basis]
    for i,j in combinations(range(14),2):
        left,right = basis[i],basis[j]
        bracket = [[sum(left[x][k]*right[k][y]-right[x][k]*left[k][y] for k in range(7))
                    for y in range(7)] for x in range(7)]
        coeff = [sign*bracket[x][y] for x,y,sign in coordinates]
        for x in range(7):
            for y in range(7):
                assert bracket[x][y] == sum(c*b[x][y] for c,b in zip(coeff,basis))
        for k,c in enumerate(coeff):
            if c:
                result[k].append((i,j,c))
    return result


def exterior_rank(degree, prime):
    target = {t:i for i,t in enumerate(combinations(range(14),degree))}
    basis = {}
    determinant = 1
    selected = []
    pivots = []
    for k,terms in enumerate(structure_constants()):
        for tail in combinations(range(14),degree-2):
            row = {}
            for i,j,c in terms:
                if i in tail or j in tail:
                    continue
                sign = (-1)**sum((i>h)+(j>h) for h in tail)
                row[target[tuple(sorted((i,j,*tail)))]] = sign*c % prime
            while row:
                lead = min(row)
                value = row[lead]
                if lead not in basis:
                    inverse = pow(value,-1,prime)
                    basis[lead] = {q:n*inverse%prime for q,n in row.items()}
                    pivots.append(lead)
                    selected.append([k,*tail])
                    determinant = determinant*value % prime
                    break
                for q,n in basis[lead].items():
                    residue = (row.get(q,0)-value*n)%prime
                    if residue:
                        row[q] = residue
                    else:
                        row.pop(q,None)
    return dict(degree=degree, prime=prime, domain=14*comb(14,degree-2),target=len(target),
                rank_lower=len(basis), determinant_in_pivot_order=determinant,
                pivot_columns=pivots, selected_rows=selected)


def partitions(total, maximum, prefix=()):
    if total==0:
        yield prefix
    for first in range(min(maximum,total),0,-1):
        yield from partitions(total-first,first,(*prefix,first))


@lru_cache(None)
def length_bound(shape, exterior_degree):
    """Find a weight in the Schur representation requiring the greatest excess."""
    best = (0, shape)
    for weight in partitions(sum(shape),shape[0]):
        if any(sum(weight[:i])>sum(shape[:i]) for i in range(1,len(weight)+1)):
            continue
        count = sum(x>exterior_degree for x in weight)+(weight.count(exterior_degree)+1)//2
        if count>best[0]:
            best = count, weight
    return best


def result():
    minor = exterior_rank(4,1000003)
    assert minor["rank_lower"]==minor["target"]==1001 and minor["determinant_in_pivot_order"]
    from g2_primary_targets import partner
    sources = [(6,5,4),(6,6,4),(6,6,6)]
    targets = []
    for source in sources:
        for order in range(1,5):
            shape = partner(source,order)
            if not shape:
                continue
            lower,weight = length_bound(shape,4)
            targets.append(dict(source=source,order=order,shape=shape,test_weight=weight,minimum_excess=lower,
                                maximum_source_length=sum(source)-1-lower))
    return dict(relation=minor,targets=targets)


if __name__ == "__main__":
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    if not __debug__:
        raise SystemExit("enable the rank checks")
    print(json.dumps(result(),separators=(",",":")))
