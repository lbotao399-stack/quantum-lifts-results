"""Pair the second residual with exact cycles containing one complementary matter fermion."""
import argparse
from functools import cache
from itertools import product
import json
from math import gcd, lcm

from colour import Q, add, compact, items, prepend, wedge
from ward import Descent

class Second:
    def __init__(self):
        self.descent = Descent()
        self.c = self.descent.c
        self.h,self.l = self.c.high[4][0],self.c.low[4,0]
        self.by_weight = {}
        for a,w in enumerate(self.c.weights):
            self.by_weight.setdefault(w,[]).append(a)

    @cache
    def eta(self,u,v,t):
        d,c = self.descent,self.c
        value = d.homogeneous(u,v,t)
        for (a,b),x in c.mu(t).items():
            uu,s = prepend(a,u)
            vv,z = prepend(b,v)
            if s and z:
                value += s*z*x*d.lambda6(uu,vv)
        return value

    @cache
    def pair(self,a,b,t):
        value = 0
        for active in ((a,b),(b,a)):
            for (r,k),x in self.c.omega(active).items():
                u,s = prepend(r,self.h)
                v,z = prepend(k,self.l)
                if s and z:
                    value += x*s*z*self.eta(u,v,t)
        return Q(value)/2

    def score(self,b,c,col):
        a,d,e = col
        value = 0
        for active,remaining in ((b,c),(c,b)):
            t,s = wedge((remaining,d,e))
            if s:
                value += s*self.pair(a,active,t)
        return value

    @cache
    def classical(self,wt):
        """Integer left relations for Q0, indexed by a matter fermion and two scalars."""
        cols = [(a,d,e) for a,d in product(range(52), repeat=2)
                for e in self.by_weight.get(tuple(k-x-y for k,x,y in zip(wt,self.c.weights[a],self.c.weights[d])),())]
        originals = []
        for a,d,e in cols:
            col = {}
            for mask,z in self.c.delta[a].items():
                p,q = items(mask)
                for i,j,sign in ((p,q,1),(q,p,-1)):
                    key = (min(d,i),max(d,i),min(e,j),max(e,j))
                    add(col,key,z*sign)
            originals.append(compact(col))
        frequency = {}
        for col in originals:
            for t in col:
                frequency[t] = frequency.get(t,0)+1
        rows = sorted(frequency,key=lambda r:(frequency[r],r))
        index = {t:i for i,t in enumerate(rows)}
        pivots,kernel = {},[]
        for i,col in enumerate(originals):
            row = {index[t]:v for t,v in col.items()}
            recipe = {i:1}
            while row:
                k = min(row)
                a = row[k]
                if k not in pivots:
                    pivots[k] = ({t:Q(x)/a for t,x in row.items()}, {t:Q(x)/a for t,x in recipe.items()})
                    break
                v,w = pivots[k]
                for t,x in v.items():
                    add(row,t,-a*x)
                    if not row[t]:
                        del row[t]
                for t,x in w.items():
                    add(recipe,t,-a*x)
                    if not recipe[t]:
                        del recipe[t]
            else:
                common = lcm(*(int(Q(x).denominator) for x in recipe.values()))
                integer = {t:int(common*x) for t,x in recipe.items()}
                divisor = gcd(*integer.values())
                integer = {t:x//divisor for t,x in integer.items()}
                rebuilt = {}
                for t,x in integer.items():
                    for j,v in originals[t].items():
                        add(rebuilt,j,x*v)
                assert not compact(rebuilt) and integer[i]
                kernel.append(integer)
        assert len(cols)==len(pivots)+len(kernel)
        return cols,len(rows),kernel

    def evaluate(self,b,c):
        weight = tuple(-x-y for x,y in zip(self.c.weights[b],self.c.weights[c]))
        cols,rows,kernel = self.classical(weight)
        # Evaluating every column also measures the size of the uncorrected test.
        scores = [self.score(b,c,col) for col in cols]
        nonzero = []
        for row in kernel:
            value = sum(x*scores[j] for j,x in row.items())
            if value:
                nonzero.append(dict(value=str(value),cycle=[[list(cols[j]),x] for j,x in sorted(row.items())]))
        return dict(first_flavour_colours=[b,c],residual_weight=weight,source_columns=len(cols),
                    classical_output_coordinates=rows,classical_rank=len(cols)-len(kernel),
                    rational_kernel=len(kernel),nonzero_raw_scores=sum(bool(x) for x in scores),
                    nonzero_on_kernel=nonzero)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('colours',type=int,nargs='+',help='consecutive pairs of first-flavour colour labels')
    args=parser.parse_args()
    if not __debug__:
        parser.error('Run with assertions enabled')
    if len(args.colours)%2 or any(not 0<=x<52 for x in args.colours):
        parser.error('Use pairs of colour labels between 0 and 51')
    search=Second()
    for b,c in zip(args.colours[::2],args.colours[1::2]):
        print(json.dumps(search.evaluate(b,c)),flush=True)
