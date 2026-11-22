"""Test predecessors of the quartic-plane one-excess source by exact bracket contractions."""
from functools import cache
from itertools import combinations, permutations
import json
from random import Random

from colour import Colour, add, compact
from six_invariants import sign


def evaluate(c, seed):
    rng = Random(seed)
    def draw():
        return {i:rng.choice((-3,-2,-1,1,2,3)) for i in rng.sample(range(52),7)}
    scalars = [draw() for _ in range(4)]
    fermions = {p:draw() for p in combinations(range(4),2)}
    def bracket(x,y):
        out = {}
        for i,a in x.items():
            for j,b in y.items():
                for k,z in c.br[i][j].items():
                    add(out,k,a*b*z)
        return compact(out)
    @cache
    def field(l, deriv):
        if len(l)==1:
            return scalars[l[0]] if not deriv else {}
        if len(l)==2:
            return bracket(scalars[l[0]],scalars[l[1]]) if deriv else fermions[l]
        if not deriv:
            return {}
        out = {}
        for j,i in enumerate(l):
            for k,v in bracket(scalars[i],fermions[l[:j]+l[j+1:]]).items():
                add(out,k,(-1)**j*v)
        return compact(out)
    def inner(x,y):
        return sum(a*b*c.metric.get((i,j),0) for i,a in x.items() for j,b in y.items())

    @cache
    def comm(x,y):
        return bracket(field(*x),field(*y))
    @cache
    def cubic(x,y,z):
        return inner(comm(x,y),field(*z))
    @cache
    def nested(x,y,z):
        return bracket(comm(x,y),field(*z))
    @cache
    def connected(x,y,z,t,u,v):
        return inner(nested(x,y,z),nested(t,u,v))
    perms = list(permutations(range(4)))
    scores = [0,0]
    for left in perms:
        for right in perms:
            # Q of the three-flavour jet, then the five undifferentiated scalars.
            group = left[:3]
            letters = [(tuple(sorted(group)), True),
                       ((left[3],), False), *(((i,), False) for i in right)]
            a,b,d,e,f,h = letters
            phase = sign(left)*sign(right)*sign(group)
            scores[0] += phase*cubic(a,b,d)*cubic(e,f,h)
            scores[1] += phase*connected(a,d,e,b,f,h)
    return scores


def calculate():
    c = Colour()
    seeds = (44330,44331)
    rows = [evaluate(c,s) for s in seeds]
    det = rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]
    assert det == -294982331904
    return dict(seeds=seeds, rank_lower=2, minor=rows, determinant=det)


if __name__=="__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
