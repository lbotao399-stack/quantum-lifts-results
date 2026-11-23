"""Check the F4 Levi strata and the colour rank of two brackets in four alternating slots."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

from colour import Colour
from source_symmetry import centralizer, character, partitions


def simple_coordinates(c):
    from flint import fmpq_mat

    simple = fmpq_mat([list(c.weights[j]) for j in c.data['raising']])
    weights = fmpq_mat([list(w) for w in c.weights]) * simple.inv()
    assert all(x.denominator == 1 for x in weights.entries())
    return [tuple(int(x) for x in row) for row in weights.tolist()]


def levi_strata(c):
    roots = simple_coordinates(c)[4:]
    positive = [r for r in roots if min(r) >= 0]
    assert len(set(positive)) == 24 and {tuple(-j for j in r) for r in positive} | set(positive) == set(roots)
    mirror = [(v, Fraction(q)) for v,q in c.data['longest_weyl']]
    assert mirror[:4] == [(i,-1) for i in range(4)]
    for i,j in combinations_with_replacement(range(52),2):
        a,ca = mirror[i]
        b,cb = mirror[j]
        left = {mirror[k][0]: v*mirror[k][1] for k,v in c.br[i][j].items()}
        right = {k:ca*cb*v for k,v in c.br[a][b].items()}
        assert left == right
    shape = (4,4,4,2)
    # Remove the vertical two-strip occupied by the antisymmetric second jet.
    pieces = []
    for small in partitions(12):
        s = small[::-1]
        if len(s) > 4:
            continue
        rest = tuple(x-y for x,y in zip(shape, s + (0,)*(4-len(s))))
        if sum(rest)==2 and all(x in (0,1) for x in rest):
            coeff = sum(Fraction(character(s,a)*character((1,1),b)*character(shape,tuple(sorted(a+b))),
                                 centralizer(a)*centralizer(b)) for a in partitions(12) for b in partitions(2))
            assert coeff == 1
            pieces.append(s)
    assert set(pieces)=={(4,3,3,2),(4,4,3,1)}
    profiles = []
    for size in range(5):
        for nodes in combinations(range(4),size):
            system = [r for r in positive if all(not n or i in nodes for i,n in enumerate(r))]
            maximum = max((sum(r) for r in system), default=0)
            central = 4-size
            profiles.append(dict(simple_nodes=nodes, positive_roots=len(system), max_height=maximum,
                                 central_rank=central,
                                 required_nilpotent_slots=[sum(mu[central:]) for mu in pieces]))
    assert max(t['max_height'] for t in profiles if t['central_rank']==0) == 11
    return dict(source_scalar_diagram=shape, scalar_coefficients=pieces,
                checked_chevalley_brackets=52*53//2, levi_subsystems=profiles)


def extend(pivots, row, p):
    row = {k:v%p for k,v in row.items() if v%p}
    while row:
        k = min(row)
        if k not in pivots:
            inverse = pow(row[k],-1,p)
            pivots[k] = {j:v*inverse%p for j,v in row.items()}
            return k
        factor = row[k]
        for j,a in pivots[k].items():
            value = (row.get(j,0)-factor*a)%p
            if value:
                row[j] = value
            else:
                row.pop(j,None)
    return None


def bracket_rank(c, p):
    pairs = list(combinations_with_replacement(range(52),2))
    order = {v:i for i,v in enumerate(pairs)}
    form = [c.metric.get((i,j),0) for i,j in pairs]
    pivots = {}
    chosen = []
    for i,j,k,l in combinations(range(52),4):
        row = Counter()
        for a,b,r,s,sign in ((i,j,k,l,1),(i,k,j,l,-1),(i,l,j,k,1)):
            for m,x in c.br[a][b].items():
                for n,y in c.br[r][s].items():
                    row[order[min(m,n),max(m,n)]] += sign*x*y
        if not row:
            continue
        assert sum(x*form[j] for j,x in row.items()) == 0
        pivot = extend(pivots,row,p)
        if pivot is not None:
            chosen.append((i,j,k,l,pivot))
    remaining = sorted(set(range(len(pairs))) - set(pivots))
    assert len(remaining)==1 and form[remaining[0]]
    from hashlib import sha256
    digest = sha256(json.dumps(chosen,separators=(',',':')).encode()).hexdigest()
    return dict(prime=p,quartets_checked=270725,symmetric_columns=len(pairs),rank=len(pivots),
                omitted_column=pairs[remaining[0]],pivot_sha256=digest)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, choices=(1009,1013), action='append')
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    c=Colour()
    result=levi_strata(c)
    result['alternating_brackets']=[bracket_rank(c,p) for p in args.prime or (1009,1013)]
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print([(r['prime'],r['symmetric_columns'],r['rank']) for r in result['alternating_brackets']])
