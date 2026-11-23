"""Lowest colour-weight-allowed transverse order at commuting F4 root quadruples."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations
import json
from pathlib import Path

from colour import Colour
from commuting_source import simple_coordinates


ZERO=(0,0,0,0)


def add_weight(a,b):
    return tuple(x+y for x,y in zip(a,b))


def reflections(c,roots):
    lookup={w:i for i,w in enumerate(roots) if i>=4}
    maps=[]
    for e,f in zip(c.data['raising'],c.data['lowering']):
        h=c.br[e][f]
        denominator=sum(c.weights[e][j]*x for j,x in h.items())
        m={}
        for j in range(4,52):
            num=2*sum(c.weights[j][i]*x for i,x in h.items())
            assert num%denominator==0
            k=num//denominator
            w=tuple(a-k*b for a,b in zip(roots[j],roots[e]))
            m[j]=lookup[w]
        assert all(m[m[j]]==j for j in m)
        maps.append(m)
    return maps


def quad_orbits(c,roots):
    positive=[i for i in range(4,52) if min(roots[i])>=0]
    allowed={v for v in combinations(positive,4) if all(not c.br[i][j] for i,j in combinations(v,2))}
    maps=reflections(c,roots)
    remain=allowed.copy()
    while remain:
        q=min(remain)
        orbit={q};todo=[q]
        while todo:
            old=todo.pop()
            for f in maps:
                new=tuple(sorted(f[j] for j in old))
                if new not in orbit:
                    orbit.add(new);todo.append(new)
        active=orbit&allowed
        assert active<=remain
        remain-=active
        yield q,len(orbit),len(active)


@cache
def occupations(k,cap):
    if len(cap)==1:
        return [(k,)] if 0<=k<=cap[0] else []
    return [(a,*v) for a in range(min(cap[0],k)+1) for v in occupations(k-a,cap[1:])]


def scalar_weights():
    from itertools import permutations
    # Permit every labelling of the four root directions and of the marked pair.
    weights=set()
    for i,j in combinations(range(4),2):
        degree=[4,4,4,2];degree[i]-=1;degree[j]-=1
        weights.update(permutations(degree))
    return sorted(weights)


class AdjointSums:
    def __init__(self,roots):
        self.step=tuple(sorted(set(roots)))
        self.levels=[{ZERO}]
    def level(self,k):
        while len(self.levels)<=k:
            self.levels.append({add_weight(v,w) for v in self.levels[-1] for w in self.step})
        return self.levels[k]
    def decomposition(self,v,k):
        if k==0:
            assert v==ZERO
            return []
        for a in self.step:
            rest=tuple(x-y for x,y in zip(v,a))
            if rest in self.levels[k-1]:
                return [a,*self.decomposition(rest,k-1)]
        raise AssertionError(v)


def calculate():
    c=Colour()
    roots=simple_coordinates(c)
    sums=AdjointSums(roots)
    weights=scalar_weights()
    records=[]
    for q,sz,positive in quad_orbits(c,roots):
        w=[roots[j] for j in q]
        for k in range(13):
            allowed=sums.level(k+1)
            hits=[]
            for degree in weights:
                for freed in occupations(k,degree):
                    remaining=tuple(sum((n-a)*w[i][j] for i,(n,a) in enumerate(zip(degree,freed))) for j in range(4))
                    if remaining in allowed:
                        hits.append((degree,freed,remaining))
            if hits:
                degree,freed,remaining=hits[0]
                witness=sums.decomposition(tuple(-x for x in remaining),k+1)
                assert all(sum(a for a in col)+b==0 for col,b in zip(zip(*witness),remaining))
                records.append(dict(root_colours=q,weyl_orbit_size=sz,positive_quartets=positive,
                                    first_allowed_order=k,degree=degree,marked_replacements=freed,
                                    compensating_adjoint_weights=witness))
                break
    by_order=defaultdict(Counter)
    for r in records:
        for t,n in (('Weyl_orbits',1),('positive_quartets',r['positive_quartets']),('all_quartets',r['weyl_orbit_size'])):
            by_order[r['first_allowed_order']][t]+=n
    assert min(by_order)>0
    return dict(weight_order='native simple roots',scalar_occupations=weights,orbit_count=len(records),
                counts_by_first_allowed_order=dict(sorted(by_order.items())),orbits=records)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',required=True,type=Path)
    a=p.parse_args()
    if not __debug__:p.error('enable assertions')
    d=calculate();a.out.write_text(json.dumps(d,indent=2)+'\n');print(d['counts_by_first_allowed_order'])
