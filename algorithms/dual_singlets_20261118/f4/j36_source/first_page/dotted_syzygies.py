"""First classical syzygies with repeated ghosts in one dotted direction and one scalar."""
import argparse
from bisect import bisect_left
from functools import cache
from math import comb
import json

from colour import add, compact, items
from mixed_syzygies import PRIME, Relations, rank_mod
from one_fermion import orbits


def insert(ghost, ordered):
    i = bisect_left(ordered,ghost)
    if i < len(ordered) and ordered[i] == ghost:
        return (),0
    return ordered[:i]+(ghost,)+ordered[i:], (-1)**i


class Dotted(Relations):
    def basis(self, degree, excess, weight):
        # One even scalar, all other letters odd. Excess counts derivatives
        # beyond the first jet on each letter. Grades here need at most two.
        ghosts = degree-excess-1
        for scalar_order in range(excess+1):
            extra = excess-scalar_order
            for twos in range(extra//2+1):
                ones = extra-2*twos
                zeros = ghosts-ones-twos
                if zeros < 0:
                    continue
                modes = (0,)*zeros+(1,)*ones+(2,)*twos
                for a,*colours in self.assignments((-1,*modes),weight):
                    tokens = tuple(52*m+c for m,c in zip(modes,colours))
                    if any(u==v for u,v in zip(tokens,tokens[1:])):
                        continue
                    yield (52*scalar_order+a,*tokens)

    @cache
    def split_scalar(self,token):
        order,colour = divmod(token,52)
        return [(52*(order-m)+j,52*(m-1)+i,sg*z)
                for m in range(1,order+1) for mask,z in self.delta[colour].items()
                for a,b in (items(mask),) for i,j,sg in ((a,b,1),(b,a,-1))]

    @cache
    def split_ghost(self,token):
        index,colour = divmod(token,52)
        order = index+1
        result = []
        for m in range(1,order//2+1):
            n = order-m
            for mask,z in self.delta[colour].items():
                a,b = items(mask)
                result.append((52*(m-1)+a,52*(n-1)+b,z))
                if m!=n:
                    result.append((52*(m-1)+b,52*(n-1)+a,-z))
        return result

    def differential(self,word):
        x,ghosts = word[0],word[1:]
        result = {}
        for y,g,z in self.split_scalar(x):
            terms,sign = insert(g,ghosts)
            if sign:
                add(result,(y,*terms),sign*z)
        for i,token in enumerate(ghosts):
            if token<52:
                continue
            spectators = ghosts[:i]+ghosts[i+1:]
            for a,b,z in self.split_ghost(token):
                second,s = insert(b,spectators)
                if s:
                    first,t = insert(a,second)
                    if t:
                        add(result,(x,*first),(-1)**i*s*t*z)
        return compact(result)

    def graviton(self,remaining):
        # 2 B(phi,c_2) + B(D phi,c_1) in divided jets.
        row = {}
        for (a,b),z in self.inverse.items():
            for x,g,v in ((a,52+b,2*z),(52+a,b,z)):
                terms,sign = insert(g,remaining)
                if sign:
                    add(row,(x,*terms),sign*v)
        return compact(row)

    def slice(self,degree,w):
        cols = list(self.basis(degree,1,w))
        index = {col:i for i,col in enumerate(cols)}
        assert len(index)==len(cols)
        images = [self.differential(col) for col in cols]
        boundaries = [self.differential(col) for col in self.basis(degree,2,w)]
        grav = []
        for colours in self.assignments((0,)*(degree-3),w):
            if len(set(colours))==len(colours):
                grav.append(self.graviton(colours))
        cycles = []
        for row in boundaries+grav:
            vector = {index[word]:c for word,c in row.items()}
            check = {}
            for i,a in vector.items():
                for word,b in images[i].items():
                    add(check,word,a*b)
            assert not compact(check)
            cycles.append(vector)
        outgoing = rank_mod(images)
        incoming = rank_mod(cycles[:len(boundaries)])
        joint = rank_mod(cycles)
        assert outgoing+joint<=len(cols)
        return dict(columns=len(cols),classical_rank_lower=outgoing,boundary_rank_lower=incoming,
                    known_cycle_rank_lower=joint,unresolved_upper=len(cols)-outgoing-joint)


def calculate(degree,weight=None):
    algebra = Dotted()
    rows = []
    for w,size in orbits(algebra.c,degree-1):
        if weight is None or w==weight:
            rows.append(dict(weight=w,orbit=size,**algebra.slice(degree,w)))
    fields = ("columns","classical_rank_lower","boundary_rank_lower","known_cycle_rank_lower","unresolved_upper")
    total = {k:sum(row["orbit"]*row[k] for row in rows) for k in fields}
    if weight is None:
        assert total["columns"] == 52*comb(52,degree-2)+52**2*comb(52,degree-3)
    else:
        assert len(rows)==1
    return dict(dotted_degree=degree-1,scalar_degree=1,prime=PRIME,blocks=sorted(rows,key=lambda v:v["weight"]),total=total)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("degree",type=int,choices=(3,4,5))
    parser.add_argument("--weight",type=int,nargs=4)
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    print(json.dumps(calculate(args.degree,tuple(args.weight) if args.weight else None),indent=2))
