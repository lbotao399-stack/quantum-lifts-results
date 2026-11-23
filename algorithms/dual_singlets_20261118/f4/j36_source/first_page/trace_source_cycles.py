"""Spin and flavour singlet trace sources modulo inactive ghost splittings at F4 level 36."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, product
import json
from math import lcm


from ghost_trace_quotient import rotate, sign
from trace_ideal import difference, echelon, monomials, relations, order, trace


PAIRS = tuple(combinations(range(3),2))
WEIGHT = (4,4,2,2,2)
GLUE = (3,4)


def odd(k):
    return k >= 3


@cache
def factor(w):
    if max(w) < 5:
        return trace(w)
    ghosts = tuple(k for k in w if k in (3,4))
    if len(set(ghosts)) != len(ghosts) or len(w)<2:
        return 0, ()
    phase = sign(ghosts)
    i = next(i for i,k in enumerate(w) if k>=5)
    phase *= (-1)**sum(k in (3,4) for k in w[i+1:])
    rest = tuple(k for k in w if k not in (3,4))
    a,b = rotate(rest),rotate(rest[::-1])
    if a==b and len(w)%2:
        return 0, ()
    return phase*(-1 if b<a and len(w)%2 else 1),tuple(sorted(ghosts))+min(a,b)


@cache
def canonical(traces):
    phase, packed = 1,[]
    for w in traces:
        s,f = factor(w)
        if not s:
            return 0, ()
        phase *= s
        packed.append(f)
    odd_factors = [f for f in packed if sum(map(odd,f))%2]
    if len(odd_factors)!=len(set(odd_factors)):
        return 0, ()
    return phase*sign(odd_factors),tuple(sorted(packed))


def sources():
    out = set()
    for i,pair in enumerate(PAIRS):
        for n in range(5):
            for scalar in product(range(3),repeat=n):
                for ghosts in ((),(3,),(4,),GLUE):
                    w = (*ghosts,5+i,*scalar)
                    s,f = factor(w)
                    if not s:
                        continue
                    wt = (ghosts.count(3),ghosts.count(4),*(scalar.count(k)+int(k in pair) for k in range(3)))
                    rest = difference(WEIGHT,wt)
                    if min(rest)<0:
                        continue
                    for tail in monomials(rest):
                        phase,word = canonical((f,*tail))
                        if phase:
                            out.add(word)
    assert all(w.count(GLUE)>=2 for w in out)
    return sorted(out)


def derivative(m,raising=None):
    out = Counter()
    before = 0
    for ix,w in enumerate(m):
        for site,k in enumerate(w):
            replacements = []
            if raising is None and k>=5:
                a,b = PAIRS[k-5]
                replacements = [((a,b),(-1)**before),((b,a),-(-1)**before)]
            elif raising:
                a,b = raising
                if k==b:
                    replacements = [((a,),1)]
                if k>=5 and a<3 and b in PAIRS[k-5] and a not in PAIRS[k-5]:
                    v = tuple(a if t==b else t for t in PAIRS[k-5])
                    replacements = [((5+PAIRS.index(tuple(sorted(v))),),sign(v))]
            for terms,c in replacements:
                f = w[:site]+terms+w[site+1:]
                s,target = canonical(m[:ix]+(f,)+m[ix+1:])
                if s:
                    out[target] += s*c
            before += odd(k)
    return {k:v for k,v in out.items() if v}


def column_matrix(columns):
    from sympy import Matrix

    keys = sorted(set().union(*(v.keys() for v in columns)))
    return Matrix(len(keys),len(columns),lambda i,j:columns[j].get(keys[i],0))


def quartic_pivots():
    all_terms = monomials(WEIGHT)
    index = {w:i for i,w in enumerate(all_terms)}
    basis = {}
    for wt,rels in relations(Fraction(1,12)).items():
        rest = difference(WEIGHT,wt)
        if min(rest)<0:
            continue
        for tail in monomials(rest):
            for rel in rels:
                row = Counter()
                for m,v in rel:
                    s,key = order(m+tail)
                    if s:
                        row[index[key]] += s*v
                echelon(basis,row)
    return basis,index


def reduce_tree(m,basis,index):
    row = {index[k]:v for k,v in derivative(m).items()}
    for j in sorted(basis):
        x = row.get(j,0)
        if x:
            for k,v in basis[j].items():
                row[k] = row.get(k,0)-x*v
    return {k:v for k,v in row.items() if v}


def coordinates():
    from sympy import Matrix, zeros

    candidates = sources()
    ward = [column_matrix([derivative(m,op) for m in candidates]) for op in ((3,4),(0,1),(1,2))]
    invariant = Matrix.hstack(*Matrix.vstack(*ward).nullspace())
    assert Matrix.vstack(*ward)*invariant == zeros(sum(m.rows for m in ward),invariant.cols)
    piv,index = quartic_pivots()
    tree = column_matrix([reduce_tree(m,piv,index) for m in candidates])*invariant
    before = column_matrix([derivative(m) for m in candidates])*invariant
    null = Matrix.hstack(*tree.nullspace())
    cycles = invariant*null
    assert tree*null == zeros(tree.rows,null.cols)
    def pack(m):
        out = []
        for j in range(m.cols):
            den = lcm(*(int(x.q) for x in m[:,j]))
            out.append([[i,int(den*x)] for i,x in enumerate(m[:,j]) if x])
        return out
    return dict(weight=WEIGHT,ghost_quotient_monomials=len(candidates),
                trace_power_distribution=dict(sorted(Counter(w.count(GLUE) for w in candidates).items())),
                singlets=invariant.cols,tree_rank_before_quartic=before.rank(),tree_rank_after_quartic=tree.rank(),
                tree_cycles=cycles.cols,monomials=candidates,singlet_columns=pack(invariant),
                cycle_columns=pack(cycles))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(coordinates(),indent=2))
