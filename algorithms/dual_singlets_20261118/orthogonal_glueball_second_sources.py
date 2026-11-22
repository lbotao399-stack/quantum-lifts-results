"""Evaluate orthogonal determinant relations on every commuting monomial and the second dual."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import factorial, prod
from pathlib import Path

from flint import fmpq_mat
from orthogonal_glueball_first_pair import ENTRIES, GHOSTS, PHI, PfaffianCoefficient
from orthogonal_glueball_second_pair import S
from terminal_commutators import words


MATRICES=ENTRIES+(((0,0,1),(1,1,1)),)
VAC={(0,0,0):1}


def odd_swaps(seq):
    return (-1)**sum(x>y for i,x in enumerate(seq) for y in seq[i+1:])


def specs(power):
    g=power-2
    patterns=[('pm','pm'),('pm','ss'),('ps','ms'),('pm','ps'),('pm','ms'),('pms',),('pm',),('closed',)]
    for pattern in patterns:
        nplus=g-sum('p' in s or s=='closed' for s in pattern)
        nminus=g-sum('m' in s or s=='closed' for s in pattern)
        if not 0<=nplus<=3 or not 0<=nminus<=3:
            continue
        num_scalar=sum(s.count('s') for s in pattern if s!='closed')
        for indices in product(range(3),repeat=num_scalar+int('pms' in pattern)):
            for plus in combinations(range(3),nplus):
                for minus in combinations(range(3),nminus):
                    yield dict(pattern=pattern,indices=indices,plus=plus,minus=minus)


class Relation:
    def __init__(self,power,spec):
        self.closed=spec['pattern']==('closed',)
        self.project=PfaffianCoefficient(2,0)
        self.catalog=defaultdict(list)
        self.mask=sum(1<<i for i in spec['plus'])+sum(1<<(3+i) for i in spec['minus'])
        self.used=0
        for i in range(3):
            for s in range(3):
                self.catalog[i,].append((s,0,1<<(3*(3*i+s)),0,1))
        for i,selected in enumerate((spec['plus'],spec['minus'])):
            for s in selected:
                self.catalog[3+i,].append((s,1<<(3*i+s),0,0,1))
        directions=iter(spec['indices'])
        for k,kind in enumerate(spec['pattern']):
            if kind=='closed':
                continue
            choices=[]
            for t in kind:
                if t=='s':
                    s=next(directions)
                    choices.append([(i,1<<(3*(3*i+s))) for i in range(3)])
                else:
                    choices.append([(3 if t=='p' else 4,0)])
            matrix=next(directions) if len(kind)==3 else 3
            bit=1<<(6+k) if sum(c!='s' for c in kind)%2 else 0
            self.mask|=bit
            self.used|=1<<k
            for fields in product(*choices):
                for order in permutations(range(len(fields))):
                    sign=odd_swaps([j for j in order if fields[j][0]>=3])
                    key=tuple(fields[j][0] for j in order)
                    degree=sum(x for _,x in fields)
                    self.catalog[key].append((matrix,bit,degree,1<<k,sign))
        self.letters=2*power+2
        self.pf_degree=self.letters-sum(len(s)-1 if s!='closed' else 2 for s in spec['pattern'])
        self.trace=lru_cache(None)(self.trace)
        self.site=lru_cache(None)(self.site)
        self.block=lru_cache(None)(self.block)
        self.raw=lru_cache(None)(self.raw)
        self.product=lru_cache(None)(self.product)
        self.commuting_log=lru_cache(None)(self.commuting_log)

    def compatible(self,a,b):
        m,x,z=a
        n,y,t=b
        if m&n or z&t or x+y not in self.project.divisors:
            return None
        return (m|n,x+y,z|t),(-1)**sum((m>>(j+1)).bit_count() for j in range(8) if n&(1<<j))

    def mul(self,a,b):
        out=Counter()
        for i,c in a.items():
            for j,d in b.items():
                match=self.compatible(i,j)
                if match:
                    k,s=match
                    out[k]+=c*d*s
        return {k:v for k,v in out.items() if v}

    def product(self,items,cartan):
        if not items:
            return VAC
        last=self.commuting_log if cartan=='ordered' else self.site if cartan else self.trace
        return self.mul(self.product(items[:-1],cartan),last(items[-1]))

    def finish(self,poly):
        sign=(-1)**(self.mask.bit_count()*(self.mask.bit_count()-1)//2)
        return sign*sum(poly.get((self.mask,x,self.used),0)*c for x,c in self.project.weights.items())

    def raw(self,word):
        """Trace-log coefficient of an ordered, uncanonicalized word."""
        states=defaultdict(Counter)
        states[0].update({(0,0,0,0,0,0):1,(1,1,0,0,0,0):1})
        for offset in range(len(word)):
            for size in (1,2,3):
                if size>len(word)-offset:
                    continue
                for b,n,y,t,c in self.catalog.get(word[offset:offset+size],()):
                    for (start,end,m,x,z,count),v in states[offset].items():
                        if not v:
                            continue
                        match=self.compatible((m,x,z),(n,y,t))
                        if match:
                            (mask,weight,used),s=match
                            for i,j,a in MATRICES[b]:
                                if i==end:
                                    states[offset+size][start,j,mask,weight,used,count+1]+=v*a*c*s
        out=Counter()
        for (i,j,m,x,z,count),v in states[len(word)].items():
            if i==j and v:
                out[m,x,z]+=Fraction((-1)**(count+1)*v,2*count)
        return out

    def trace(self,word):
        from mixed_incoming import so_canon
        alphabet=PHI+GHOSTS
        out=Counter()
        for raw in {v[i:]+v[:i] for v in (word,word[::-1]) for i in range(len(word))}:
            s,cw=so_canon((tuple(alphabet[j] for j in raw),))
            assert cw==(tuple(alphabet[j] for j in word),) and s
            for k,v in self.raw(raw).items():
                out[k]+=s*v
        return {k:v for k,v in out.items() if v}

    def coefficient(self,w):
        if self.closed:
            pure=tuple((PHI+GHOSTS).index(l) for l in S)
            if pure not in w:
                return Fraction()
            w=list(w)
            w.remove(pure)
            w=tuple(w)
        return Fraction(self.finish(self.product(w,False)),prod(factorial(n) for n in Counter(w).values()))

    def block(self,word):
        # Convert every ordered colour path to a graded-commuting Cartan
        # coefficient with a single standard 2-by-2 colour rotation.
        out=Counter()
        for raw in set(permutations(word)):
            s=odd_swaps([j for j in raw if j>=3])
            for b,m,x,z,c in self.catalog.get(raw,()):
                out[b,m,x,z]+=s*c*(-1)**(len(word)//2)
        return {k:v for k,v in out.items() if v}

    def site(self,weight):
        word=tuple(j for j,k in enumerate(weight[2:]+weight[:2]) for _ in range(k))
        if not word:
            return VAC
        out=Counter()
        for (b,m,x,z),c in self.block(word).items():
            if b==3:
                out[m,x,z]+=2*c
        for size in range(1,len(word)):
            if size>3 or len(word)-size>3:
                continue
            for left in set(combinations(word,size)):
                right=list(word)
                for j in left:
                    right.remove(j)
                right=tuple(right)
                orient=odd_swaps([j for j in left+right if j>=3])
                for (b,m,x,z),c in self.block(left).items():
                    for (a,n,y,t),d in self.block(right).items():
                        if (a==3)!=(b==3):
                            continue
                        metric=sum(v*w for i,j,v in MATRICES[a] for l,k,w in MATRICES[b] if i==k and j==l)
                        match=self.compatible((m,x,z),(n,y,t))
                        if metric and match:
                            key,s=match
                            out[key]+=Fraction(c*d*metric*s*orient,2)
        return {k:v for k,v in out.items() if v}

    def cartan(self,sites):
        if not self.closed:
            return self.finish(self.product(sites,True))
        value=0
        for k,w in enumerate(sites):
            if w[:2]==(1,1):
                modified=sites[:k]+((0,0)+w[2:],)+sites[k+1:]
                value+=2*self.finish(self.product(modified,True))
        return value

    def commuting_log(self,weight):
        from orthogonal_glueball_first_pair import multiset_orders
        word=tuple(j for j,k in enumerate(weight[2:]+weight[:2]) for _ in range(k))
        out=Counter()
        for raw in multiset_orders(word):
            sign=odd_swaps([x for x in raw if x>=3])
            for key,c in self.raw(raw).items():
                out[key]+=sign*c
        return {k:v for k,v in out.items() if v}

    def commuting_coefficient(self,traces):
        sign=1
        if self.closed:
            pure=(1,1,0,0,0)
            if pure not in traces:
                return Fraction()
            work=list(traces)
            work.remove(pure)
            traces=tuple(work)
            sign=-1
        return Fraction(sign*self.finish(self.product(traces,'ordered')),
                        prod(factorial(n) for n in Counter(traces).values()))

    def clear(self):
        for method in (self.raw,self.product,self.block,self.site,self.trace,self.commuting_log):
            method.cache_clear()


def evaluate(args):
    spec,dual,colours,power=args
    relation=Relation(power,spec)
    assert relation.pf_degree>colours
    cartan=[str(relation.cartan(v)) for v in words((power-2,power-2,2,2,2))]
    score=sum(Fraction(c)*relation.coefficient(tuple(tuple(t) for t in w)) for w,c in dual)
    relation.clear()
    return dict(spec=spec,pfaffian_degree=relation.pf_degree,cartan=cartan,second_score=str(score))


def ranks(rows):
    matrix=fmpq_mat([row['cartan'] for row in rows])
    augmented=fmpq_mat([row['cartan']+[row['second_score']] for row in rows])
    first,second=matrix.rank(),augmented.rank()
    result=dict(cartan_rank=first,augmented_rank=second)
    if second>first:
        reduced,_=matrix.transpose().rref()
        pivots=[next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(first)]
        for j in range(len(rows)):
            if j in pivots:
                continue
            vec=[Fraction(0)]*len(rows)
            vec[j]=1
            for i,k in enumerate(pivots):
                vec[k]=-Fraction(str(reduced[i,j]))
            assert not any(sum(c*Fraction(rows[z]['cartan'][k]) for z,c in enumerate(vec))
                           for k in range(matrix.ncols()))
            score=sum(c*Fraction(rows[z]['second_score']) for z,c in enumerate(vec))
            if score:
                result['witness']={str(z):str(c) for z,c in enumerate(vec) if c}
                result['witness_score']=str(score)
                break
        assert 'witness' in result
    return result


if __name__=='__main__':
    from concurrent.futures import ProcessPoolExecutor
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dual',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--range',nargs=2,type=int)
    a=p.parse_args()
    if not __debug__:
        p.error('Run with assertions enabled')
    d=json.loads(a.dual.read_bytes())
    n,power=d['colours'],d['power']
    tasks=[(s,d['scores'],n,power) for s in specs(power) if Relation(power,s).pf_degree>n]
    if a.range:
        tasks=tasks[slice(*a.range)]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows=list(pool.map(evaluate,tasks))
    result=dict(colours=n,power=power,rows=rows,**ranks(rows))
    a.out.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(len(rows),result['cartan_rank'],result['augmented_rank'],result.get('witness_score'))
