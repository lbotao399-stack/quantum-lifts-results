"""Evaluate the dual descent from native brackets on all extremal scalar weights."""
import argparse
from collections import Counter
from functools import cache, lru_cache
from itertools import combinations
import json
from math import prod

from colour import Colour, Q, add, apply, compact, items, prepend, wedge

class Descent:
    def __init__(self):
        self.c = c = Colour()
        self.high_weight = {n:[c.weight(h) for h in group] for n,group in c.high.items()}
        self.low_weight = {n:[c.weight(c.low[n,i]) for i in range(len(group))] for n,group in c.high.items()}
        self.mirror = [(i,Q(v)) for i,v in c.data['longest_weyl']]
        self.phase = [prod(self.mirror[i][1] for i in items(h)) for h in c.high[5]]
        self.pinned = {}
        self.seeds = []
        coefficients = iter(map(Q,c.data['coefficients']))
        for section in c.data['seeds']:
            i,j = section['j'],section['i']
            out = {}
            for entries in section['basis']:
                vector = {wedge((a,b,k))[0]:v for a,b,k,v in entries}
                assert all(a<b<k for a,b,k,v in entries)
                self.seeds.append((i,j,vector))
                scalar = next(coefficients)
                for m,v in vector.items():
                    add(out,m,scalar*v)
            self.pinned[i,j] = compact(out)
        assert next(coefficients,None) is None and len(self.seeds)==40
        self.colour_duals = [[(i,c.metric[i,j]) for i in range(52) if (i,j) in c.metric] for j in range(52)]
        self.lambda_weights = self.classical_weights()

    @cache
    def dual(self,t):
        row = {0:1}
        for i in reversed(items(t)):
            new = {}
            for v,x in row.items():
                for j,y in self.colour_duals[i]:
                    w,s = prepend(j,v)
                    add(new,w,s*x*y)
            row = compact(new)
        return row

    @cache
    def denominator(self,n,w,raising):
        shift = self.c.rho2
        other = tuple(a+(b if raising else -b) for a,b in zip(w,shift))
        den = 18*n - self.c.scalar_product(w,other)
        return -Q(2)/den if den>0 else 0

    @lru_cache(260000)
    def pure(self,i,u,v):
        c = self.c
        weight = c.weight(u)
        if weight in self.high_weight[6]:
            if u == c.high[6][i] and v == c.low[6,i]:
                return c.norm[6,i]*c.low_sign[6,i]
            return 0
        factor = self.denominator(6,weight,True)
        if not factor:
            return 0
        value = 0
        for e,f,inv in c.roots:
            raised = c.act(e,u)
            if raised:
                for vv,y in c.act(f,v).items():
                    for uu,x in raised.items():
                        value += inv*x*y*self.pure(i,uu,vv)
        return value*factor

    def classical_weights(self):
        # A mixed bracket on the two old five-ghost legs factors through their
        # abelian-ideal modules. Test every weight of the remaining adjoint leg.
        rows = []
        for j,h in enumerate(self.c.high[5]):
            for k in range(2):
                l = self.c.low[5,k]
                for r,w in enumerate(self.c.weights):
                    if any(a+b+d for a,b,d in zip(self.c.weight(h),self.c.weight(l),w)):
                        continue
                    vector = [0,0]
                    for pair,z in self.c.delta[r].items():
                        a,b = items(pair)
                        for x,y,coef in ((a,b,z),(b,a,-z)):
                            u,s = prepend(x,h)
                            v,t = prepend(y,l)
                            if s and t:
                                for i in range(2):
                                    vector[i] += s*t*coef*self.pure(i,u,v)
                    rows.append((j,k,r,vector))
        self.classical_rows = rows
        x,y = next(row for j,k,r,row in rows if any(row))
        from math import gcd,lcm
        d = lcm(int(Q(x).denominator),int(Q(y).denominator))
        p,q = int(d*y),int(-d*x)
        g = gcd(p,q) * (1 if p>0 else -1)
        weights = p//g,q//g
        assert weights == (2392,2299) and all(sum(a*b for a,b in zip(row,weights))==0 for *_,row in rows)
        return weights

    @cache
    def lambda6(self,u,v):
        c = self.c
        if c.scalar_product(c.weight(u),c.rho2) < 0:
            u,v = v,u
        return sum(k*self.pure(i,u,v) for i,k in enumerate(self.lambda_weights))

    @cache
    def transformed(self,t):
        terms = [self.mirror[i] for i in items(t)]
        u,s = wedge([i for i,_ in terms])
        return u,s*prod(x for _,x in terms)

    @cache
    def pin(self,i,j,t):
        if (i,j)==(0,1):
            u,z = self.transformed(t)
            return z*self.phase[0]/self.phase[1]*self.pin(1,0,u)
        return sum(v*self.pinned[i,j].get(u,0) for u,v in self.dual(t).items())

    @lru_cache(320000)
    def lower(self,i,j,v,t):
        c = self.c
        w = c.weight(v)
        if any(a+b+d for a,b,d in zip(w,c.weight(t),self.high_weight[5][i])):
            return 0
        if w in self.low_weight[5]:
            return c.low_sign[5,j]*self.pin(i,j,t) if v==c.low[5,j] else 0
        factor = self.denominator(5,w,False)
        if not factor:
            return 0
        value = 0
        for e,f,g in c.roots:
            lowered = c.act(f,v)
            if lowered:
                for tt,y in c.act(e,t).items():
                    for vv,x in lowered.items():
                        value += g*x*y*self.lower(i,j,vv,tt)
        return factor*value

    @lru_cache(400000)
    def upper(self,i,j,u,v,t):
        c = self.c
        w = c.weight(u)
        if w in self.high_weight[5]:
            return self.lower(i,j,v,t) if u==c.high[5][i] else 0
        factor = self.denominator(5,w,True)
        if not factor:
            return 0
        value = 0
        for e,f,g in c.roots:
            raised = c.act(e,u)
            if raised:
                for uu,x in raised.items():
                    for vv,y in c.act(f,v).items():
                        value += g*x*y*self.upper(i,j,uu,vv,t)
                    for tt,y in c.act(f,t).items():
                        value += g*x*y*self.upper(i,j,uu,v,tt)
        return factor*value

    def homogeneous(self,u,v,t):
        c = self.c
        if c.scalar_product(c.weight(v),c.rho2)>c.scalar_product(c.weight(u),c.rho2):
            u,v = v,u
        return sum(self.upper(i,j,u,v,t) for i in range(2) for j in range(2))

    def dphi(self,i,l,r,a,b):
        c = self.c
        h = c.high[5][i]
        active = wedge((a,b))[0]
        particular = homogeneous = 0
        for pair,z in c.delta[r].items():
            x,y = items(pair)
            for first,second,coef in ((x,y,z),(y,x,-z)):
                v,sv = prepend(first,l)
                t,st = prepend(second,active)
                if not (sv and st):
                    continue
                coef *= sv*st
                homogeneous += coef*sum(self.lower(i,j,v,t) for j in range(2))
                for (s,k),val in c.mu(t).items():
                    uu,su = prepend(s,h)
                    vv,ss = prepend(k,v)
                    if su and ss:
                        particular += coef*val*su*ss*self.lambda6(uu,vv)
        return particular,homogeneous

    def mixed(self,r,t):
        c = self.c
        h,l = c.high[4][0],c.low[4,0]
        result = 0
        for pair,z in c.delta[r].items():
            a,b = items(pair)
            for first,second,coef in ((a,b,z),(b,a,-z)):
                u,s = prepend(first,h)
                v,ss = prepend(second,l)
                if s and ss:
                    result += coef*s*ss*self.homogeneous(u,v,t)
        return result

    def check_seeds(self):
        c = self.c
        # The longest Weyl lift must preserve the native bracket, including
        # the phases used to exchange the two five-ghost slots.
        for a in range(52):
            b,phase = self.mirror[a]
            assert c.weights[b] == tuple(-x for x in c.weights[a])
            for j in range(52):
                k,z = self.mirror[j]
                lhs = {self.mirror[t][0]:v*self.mirror[t][1] for t,v in c.br[a][j].items()}
                rhs = {t:phase*z*v for t,v in c.br[b][k].items()}
                assert lhs==rhs and phase*z*c.metric.get((b,k),0)==c.metric.get((a,j),0)
        labels = ((0,1,0,3),(0,0,3,0))
        for i,j,vector in self.seeds:
            expect = tuple(a-b for a,b in zip(self.high_weight[5][i],self.high_weight[5][j]))
            assert all(c.weight(v)==expect for v in vector)
            for roots,powers in ((c.data['raising'],labels[j]),(c.data['lowering'],labels[i])):
                for root,power in zip(roots,powers):
                    row = vector
                    for k in range(power+1):
                        row = apply(lambda t:c.act(root,t),row)
                    assert not row
            if i==j:
                reflected = {}
                for t,x in vector.items():
                    u,z = self.transformed(t)
                    add(reflected,u,z*x)
                assert reflected==vector
        return len(self.seeds)

    def target_value(self):
        c = self.c
        dimensions = []
        for h in c.high[6]:
            weight = c.weight(h)
            shifted = tuple(2*a+b for a,b in zip(weight,c.rho2))
            n = prod(c.scalar_product(shifted,c.weights[r])/c.scalar_product(c.rho2,c.weights[r])
                     for r in c.positive)
            assert n.denominator==1
            dimensions.append(int(n))
        from math import factorial
        target = -factorial(6)*sum(a*b for a,b in zip(dimensions,self.lambda_weights))
        assert target == -15313094596800
        return dimensions,target

    def check_psi(self):
        c = self.c
        assert c.scalar_block() == 1378
        for i in range(52):
            for j in range(52):
                out = {}
                for mask,x in c.delta[i].items():
                    t,s = prepend(j,mask)
                    if s:
                        add(out,t,x*s)
                got = apply(c.mu,compact(out))
                # Matrix element of the native symmetric one-loop colour fan.
                expected = {}
                for a,b in ((i,j),(j,i)):
                    for tup,z in c.omega((a,b)).items():
                        add(expected,tup,z/2)
                assert got == compact(expected), (i,j)
        for i,j,seed in self.seeds:
            contraction = Counter()
            for t,z in seed.items():
                a,b,d = items(t)
                for x,y,s in ((a,b,d),(b,d,a),(d,a,b)):
                    for k,coef in c.br[x][y].items():
                        contraction[k,s] += z*coef
            assert not compact(contraction)
        return dict(ordered_scalar_inputs=2704, symmetric_casimir_columns=1378, homogeneous_tensors=40)

    def check_dphi(self):
        c = self.c
        l = c.low[4,0]
        by_weight = {}
        for a,b in combinations(range(52),2):
            t = wedge((a,b))[0]
            by_weight.setdefault(c.weight(t),[]).append((a,b))
        out = []
        for i,h in enumerate(c.high[5]):
            tested = counter = 0
            for r,w in enumerate(c.weights):
                wanted = tuple(-x-y-z for x,y,z in zip(c.weight(h),c.weight(l),w))
                for a,b in by_weight.get(wanted,[]):
                    first,second = self.dphi(i,l,r,a,b)
                    assert first+second==0, (i,r,a,b,str(first),str(second))
                    counter += bool(first)
                    tested += 1
            out.append(dict(raw_weights=tested,nonzero_cancellations=counter))
        assert [r['raw_weights'] for r in out] == [326,437]
        return out

    def check_mixed(self,residues,modulus):
        by_weight = {}
        for abc in combinations(range(52),3):
            t,s = wedge(abc)
            by_weight.setdefault(self.c.weight(t),[]).append(t)
        count = selected = 0
        for r,w in enumerate(self.c.weights):
            for t in by_weight[tuple(-x for x in w)]:
                if count % modulus in residues:
                    z = self.mixed(r,t)
                    assert not z,(r,items(t),str(z))
                    selected += 1
                count += 1
        assert count == 6920
        return dict(full_weights=count, tested=selected, residues=sorted(residues),modulus=modulus)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sector',choices=('target','seeds','psi','dphi','mixed'))
    parser.add_argument('--modulus',type=int,default=1)
    parser.add_argument('--residues',type=int,nargs='+',default=[0])
    args = parser.parse_args()
    if not __debug__:
        parser.error('Run with Python assertions enabled')
    assert args.modulus>0 and all(0<=r<args.modulus for r in args.residues)
    d = Descent()
    if args.sector=='target':
        dimensions,value = d.target_value()
        out = dict(rows=[[i,j,r,list(map(str,z))] for i,j,r,z in d.classical_rows],weights=d.lambda_weights,
                   branch_dimensions=dimensions,lambda_S6=value)
    elif args.sector=='seeds':
        out = dict(verified_extremal_tensors=d.check_seeds())
    elif args.sector=='psi':
        out = d.check_psi()
    elif args.sector=='dphi':
        out = d.check_dphi()
    else:
        out = d.check_mixed(set(args.residues),args.modulus)
    print(json.dumps(out,indent=2))
