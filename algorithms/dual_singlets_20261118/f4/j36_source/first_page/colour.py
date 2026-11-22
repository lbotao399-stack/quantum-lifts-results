"""Exterior powers, extremal vectors and the one-loop scalar intertwiner in F4."""
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import prod
import gzip
import json
from pathlib import Path

try:
    from flint import fmpq as Q
except ImportError:
    from fractions import Fraction as Q

def wedge(indices):
    mask = 0
    sign = 1
    for i in indices:
        if (mask >> i) & 1:
            return 0, 0
        if (mask >> (i+1)).bit_count() & 1:
            sign = -sign
        mask |= 1 << i
    return mask, sign

@cache
def items(mask):
    if not mask:
        return ()
    bit = mask & -mask
    return (bit.bit_length()-1,) + items(mask ^ bit)

def prepend(i, mask):
    if mask & (1 << i):
        return mask, 0
    return mask | (1 << i), (-1)**((mask & ((1 << i)-1)).bit_count() & 1)

def add(row, key, value):
    if value:
        row[key] = row.get(key, 0) + value

def compact(row):
    return {k:v for k,v in row.items() if v}

def apply(operator, row):
    out = {}
    for key, x in row.items():
        for target, y in operator(key).items():
            add(out,target,x*y)
    return compact(out)

class Colour:
    def __init__(self):
        self.data = json.loads(gzip.decompress(Path(__file__).with_name('inputs.json.gz').read_bytes()))
        self.br = [[{} for j in range(52)] for i in range(52)]
        for i,j,row in self.data['brackets']:
            self.br[i][j] = dict(row)
            self.br[j][i] = {k:-v for k,v in row}
        self.metric = {(i,j):int(v) for i,j,v in self.data['metric']}
        self.inverse = {(i,j):Q(v) for i,j,v in self.data['inverse']}
        for i,j in combinations_with_replacement(range(52),2):
            assert sum(self.metric.get((i,k),0)*self.inverse.get((k,j),0) for k in range(52)) == (i==j)
        self.weights = [tuple(self.br[h][i].get(i,0) for h in range(4)) for i in range(52)]
        roots = {w:i for i,w in enumerate(self.weights) if i>=4}
        self.opposite = {i:roots[tuple(-x for x in w)] for w,i in roots.items()}
        self.positive = [i for w,i in roots.items() if w>(0,0,0,0)]
        self.roots = [(i,self.opposite[i],self.inverse[i,self.opposite[i]]) for i in self.positive]
        self.rho2 = tuple(sum(self.weights[i][j] for i in self.positive) for j in range(4))
        self.high = {}
        self.low = {}
        self.low_sign = {}
        self.norm = {}
        for degree in (4,5,6):
            higher = {i:{k for e in self.data['raising'] for k in self.br[e][i]} for i in self.positive}
            ideals = []
            for u in combinations(sorted(self.positive),degree):
                support = set(u)
                if any(higher[i] - support for i in u) or any(self.br[i][j] for i,j in combinations(u,2)):
                    continue
                h = wedge(u)[0]
                labels = []
                for e,f in zip(self.data['raising'],self.data['lowering']):
                    coroot = self.br[e][f]
                    den = sum(x*self.weights[e][k] for k,x in coroot.items())
                    labels.append(Q(2*sum(x*self.weight(h)[k] for k,x in coroot.items()),den))
                ideals.append((tuple(labels),h))
            # The archive indexes its A5 seed tensors by (0,1,0,3) first.
            # Its two A6 components use lexicographic Dynkin-label order.
            ordered = sorted(ideals,reverse=(degree==5))
            if degree==5:
                assert [lab for lab,h in ordered]==[(0,1,0,3),(0,0,3,0)]
            self.high[degree] = [h for lab,h in ordered]
            assert len(ideals) == {4:1,5:2,6:2}[degree]
            for j,h in enumerate(self.high[degree]):
                u = items(h)
                l = tuple(self.opposite[i] for i in u)
                self.low[degree,j], self.low_sign[degree,j] = wedge(l)
                self.norm[degree,j] = prod(self.metric[i,k] for i,k in zip(u,l))
                assert self.scalar_product(self.weight(h),tuple(a+b for a,b in zip(self.weight(h),self.rho2))) == 18*degree
        self.delta = [dict() for i in range(52)]
        duals = defaultdict(list)
        for (i,j),x in self.inverse.items():
            duals[i].append((j,x))
        for i,j in combinations(range(52),2):
            for k,v in self.br[i][j].items():
                for (h,t),b in self.metric.items():
                    if t!=k:
                        continue
                    for a,x in duals[i]:
                        for c,y in duals[j]:
                            key,s = wedge((a,c))
                            add(self.delta[h],key,v*b*x*y*s)
        self.delta = [compact(row) for row in self.delta]

    @cache
    def weight(self,mask):
        return tuple(sum(self.weights[i][j] for i in items(mask)) for j in range(4))

    @cache
    def scalar_product(self,u,v):
        return sum(x*y*self.inverse.get((i,j),0) for i,x in enumerate(u) for j,y in enumerate(v))

    @cache
    def act(self,root,mask):
        result = {}
        for position,i in enumerate(items(mask)):
            rest = mask ^ (1<<i)
            for j,x in self.br[root][i].items():
                dest, sign = prepend(j,rest)
                add(result,dest,(-1)**(position & 1)*sign*x)
        return compact(result)

    @cache
    def omega(self,pair):
        u,v = pair
        out = {}
        for (i,j),x in self.inverse.items():
            for a,y in self.br[i][u].items():
                for b,z in self.br[j][v].items():
                    add(out,(a,b),x*y*z)
        return compact(out)

    @cache
    def sym_casimir(self,pair):
        out = {pair:36}
        for (i,j),x in self.omega(pair).items():
            add(out,tuple(sorted((i,j))),2*x)
        return compact(out)

    @cache
    def split(self,pair):
        a,b = pair
        out = {}
        for i,j in ((a,b),(b,a)):
            for pair_mask,x in self.delta[i].items():
                mask,sign = prepend(j,pair_mask)  # a two-form commutes with j
                add(out,mask,sign*x/2)
        return compact(out)

    @cache
    def merge(self,mask):
        a,b,c = items(mask)
        out = {}
        for i,j,k in ((a,b,c),(b,c,a),(c,a,b)):
            for t,x in self.br[i][j].items():
                add(out,tuple(sorted((t,k))),x)
        return compact(out)

    def shifted(self,row,c):
        out = apply(self.sym_casimir,row)
        for key,x in row.items():
            add(out,key,-c*x)
        return compact(out)

    @cache
    def mu(self,triple):
        """Symmetric ordered colour tensor on the two emitted dotted ghosts."""
        row = self.merge(triple)
        # Interpolate the inverse of (C-54)/2 at 0,26,40.
        scalar = {}
        for c in (0,26,40):
            out,den = row,1
            for k in (0,26,40):
                if k!=c:
                    out = self.shifted(out,k)
                    den *= c-k
            for key,x in out.items():
                add(scalar,key,Q(2,den*(c-54))*x)
        out = {}
        for pair,x in scalar.items():
            for (i,j),y in self.omega(pair).items():
                add(out,(i,j),x*y/2)
                add(out,(j,i),x*y/2)
        return compact(out)

    def scalar_block(self):
        for a,b in combinations_with_replacement(range(52),2):
            row = {(a,b):1}
            merged = apply(self.merge,self.split((a,b)))
            assert {k:2*v for k,v in merged.items()} == self.shifted(row,54)
            assert not self.shifted(self.shifted(self.shifted(row,0),26),40)
        return 1378
