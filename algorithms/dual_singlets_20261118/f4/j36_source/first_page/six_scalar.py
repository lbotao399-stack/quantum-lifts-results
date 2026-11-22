"""One matter fermion and four scalars in the SL(3) singlet at degree (2,2,2)."""
import argparse
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations
import json

from colour import add, compact, items
from mixed_syzygies import PRIME, Relations, rank_mod
from one_fermion import orbits, weights_sum


def alternating(t):
    if len(set(t))!=len(t):
        return (),0
    sign = (-1)**sum(a>b for i,a in enumerate(t) for b in t[i+1:])
    return tuple(sorted(t)),sign


def minors_expanded(u,v):
    out = {}
    for a in permutations(u):
        for b in permutations(v):
            key=tuple(tuple(sorted(pair)) for pair in zip(a,b))
            add(out,key,alternating(a)[1]*alternating(b)[1])
    return compact(out)


@cache
def straightening(u,v):
    """Straighten the two 3x3 minors once for each relative order of six indices."""
    letters = tuple(sorted(u+v))
    basis = set()
    for choice in combinations(range(6),3):
        a=tuple(letters[i] for i in choice)
        b=tuple(letters[i] for i in range(6) if i not in choice)
        if len(set(a))==len(set(b))==3 and all(i<=j for i,j in zip(a,b)):
            basis.add((a,b))
    pivots = {}
    def reduce_column(row,recipe):
        while row and (k:=min(row)) in pivots:
            column,coeffs=pivots[k]
            x=row[k]
            for t,y in column.items():
                add(row,t,-x*y)
            for t,y in coeffs.items():
                add(recipe,t,-x*y)
            row=compact(row)
        return row,compact(recipe)
    for pair in sorted(basis):
        row,recipe=reduce_column(minors_expanded(*pair),{pair:1})
        assert row
        k=min(row)
        x=row[k]
        pivots[k]=({t:Fraction(y,x) for t,y in row.items()},{t:Fraction(y,x) for t,y in recipe.items()})
    original=minors_expanded(u,v)
    row,recipe=reduce_column(original.copy(),{})
    assert not row
    terms=[]
    check={}
    for (a,b),y in recipe.items():
        assert y.denominator==1
        z=-int(y)
        terms.append((a+b,z))
        for t,x in minors_expanded(a,b).items():
            add(check,t,z*x)
    assert compact(check)==original
    return tuple(terms)


def packed(t):
    a,b,c,d,e,f=t
    return (((((a*64+b)*64+c)*64+d)*64+e)*64+f)


def product_rows(u,v):
    if all(a<=b for a,b in zip(u,v)):
        return ((packed(u+v),1),)
    if all(a>=b for a,b in zip(u,v)):
        return ((packed(v+u),1),)
    labels=sorted(set(u+v))
    order={a:i for i,a in enumerate(labels)}
    a,b=(tuple(order[k] for k in col) for col in (u,v))
    return [(packed(tuple(labels[i] for i in indices)),z) for indices,z in straightening(a,b)]


def source_terms(a,b,triple):
    """Standard tableaux of column lengths three and one, with a separate fermion colour."""
    col,sg=alternating(triple)
    if not sg:
        return ()
    c,d,e=col
    if b>=c:
        return (((a,b,c,d,e),sg),)
    return (((a,c,b,d,e),sg),((a,d,b,c,e),-sg),((a,e,b,c,d),sg))


class Six(Relations):
    def __init__(self):
        super().__init__()
        self.pairs=defaultdict(list)
        self.triples=defaultdict(list)
        for a in range(52):
            for b in range(52):
                self.pairs[weights_sum(self.c.weights[a],self.c.weights[b])].append((a,b))
        for triple in combinations(range(52),3):
            self.triples[weights_sum(*(self.c.weights[i] for i in triple))].append(triple)

    @cache
    def tree_minor(self,a,b):
        row={}
        for mask,z in self.delta[a].items():
            tr,sg=alternating((b,*items(mask)))
            if sg:
                add(row,tr,sg*z)
        return compact(row)

    def block(self,w):
        def complement(a):
            return tuple(t-i for t,i in zip(w,a))
        cols=[(a,b,*t) for weight,ps in self.pairs.items() for a,b in ps
              for t in self.triples.get(complement(weight),()) if b>=t[0]]
        idx={c:i for i,c in enumerate(cols)}
        assert len(idx)==len(cols)
        images=[]
        for a,b,c,d,e in cols:
            row={}
            for tr,z in self.tree_minor(a,b).items():
                for k,y in product_rows(tr,(c,d,e)):
                    add(row,k,z*y)
            images.append(compact(row))

        def append(row,a,b,t,value):
            for k,sg in source_terms(a,b,t):
                add(row,idx[k],sg*value)
        boundaries=[]
        for weight,ps in self.pairs.items():
            other=self.pairs.get(complement(weight),())
            for a,b in ps:
                for c,d in other:
                    if (a,b)>=(c,d):
                        continue
                    row={}
                    for f,g,h,i,sg in ((a,b,c,d,1),(c,d,a,b,-1)):
                        for t,z in self.tree_minor(f,g).items():
                            append(row,h,i,t,sg*z)
                    boundaries.append(compact(row))
        for a,wa in enumerate(self.c.weights):
            for tr in self.triples.get(complement(wa),()):
                row={}
                for mask,z in self.delta[a].items():
                    i,j=items(mask)
                    append(row,i,j,tr,z)
                    append(row,j,i,tr,-z)
                boundaries.append(compact(row))
        products=[]
        for b,wb in enumerate(self.c.weights):
            for c,d in self.pairs.get(complement(wb),()):
                if c>=d:
                    continue
                row={}
                for (a,i),z in self.inverse.items():
                    append(row,a,b,(i,c,d),3*z)
                    append(row,a,i,(b,c,d),-z)
                products.append(compact(row))
        for v in boundaries+products:
            check={}
            for j,x in v.items():
                for key,y in images[j].items():
                    add(check,key,x*y)
            assert not compact(check)
        tree=rank_mod(images)
        boundary=rank_mod(boundaries)
        joint=rank_mod(boundaries+products)
        assert tree+joint<=len(cols)
        return dict(columns=len(cols),classical_rank_lower=tree,boundary_rank_lower=boundary,
                    known_cycle_rank_lower=joint,unresolved_upper=len(cols)-tree-joint)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--weight",type=int,nargs=4)
    args=p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    engine=Six()
    weights=orbits(engine.c,5)
    blocks=[dict(weight=w,orbit=n,**engine.block(w)) for w,n in weights
            if args.weight is None or tuple(args.weight)==w]
    fields=("columns","classical_rank_lower","boundary_rank_lower","known_cycle_rank_lower","unresolved_upper")
    total={k:sum(r["orbit"]*r[k] for r in blocks) for k in fields}
    if args.weight is None:
        assert total["columns"]==52*(52*53*51*50//8)
    else:
        assert len(blocks)==1
    print(json.dumps(dict(flavour_degrees=[2,2,2],flavour_singlet=True,prime=PRIME,
                          blocks=sorted(blocks,key=lambda r:r["weight"]),total=total),indent=2))
