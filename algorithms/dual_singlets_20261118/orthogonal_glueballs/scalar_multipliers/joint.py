"""Verify three transported families on the remaining SO(13) second coordinate."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import json
from itertools import combinations, combinations_with_replacement, product

from flint import fmpq_mat
from second import test
from verify import GHOSTS, HERE, PHI, determinant, encode, path_types, scalar_projector, so_canon, unit, words
from terminal_commutators import parity


def mixed_specs():
    weight = (0,1,-1,0)
    for paths,p,m,_ in path_types(4,4,6,2):
        closed = [v for c,v in paths if c]
        opens = [v for c,v in paths if not c]
        assert closed in ([],['mp'])
        for matrices in product(*(range(3) if len(v)%2 else [3] for v in opens)):
            for scalar in product(range(3),repeat=sum(v.count('s') for v in opens)):
                slots = iter(scalar)
                directions = [[unit(3,next(slots)) for _ in range(v.count('s'))] for v in opens]
                for ps,ms in product(combinations(range(3),p),combinations(range(3),m)):
                    for axes in combinations_with_replacement(range(3),2):
                        charge = sum(weight[i] for i in ps+ms+axes+matrices) - sum(weight[i] for i in scalar)
                        if not charge:
                            yield dict(closed=closed,paths=opens,matrices=[unit(4,i) for i in matrices],
                                       single=[list(ps),list(ms)],scalar_directions=directions),list(axes)


def ghost_form(row, scalar):
    plus,minus = GHOSTS
    paths = ((plus,scalar,minus),(minus,scalar,plus))
    out = Counter()
    for word,c in row.items():
        candidates = set()
        for i,tr in enumerate(word):
            if len(tr) < 4:
                continue
            for j in range(len(tr)):
                rotation = tr[j:]+tr[:j]
                if rotation[:3] in paths:
                    sign,v = so_canon(word[:i]+((scalar,)+rotation[3:],)+word[i+1:])
                    if sign:
                        candidates.add(v)
        for v in candidates:
            for i,tr in enumerate(v):
                for j,old in enumerate(tr):
                    if old != scalar:
                        continue
                    for sign,seq in zip((1,-1),paths):
                        phase,w = so_canon(v[:i]+(tr[:j]+seq+tr[j+1:],)+v[i+1:])
                        if phase and w == word:
                            out[v] += sign*phase*c
    return [(encode(w),c) for w,c in out.items() if c]


def base_relation(record):
    spec = record['spec']
    rel = determinant(6-len(spec['closed']),spec)
    rel.closed = bool(spec['closed'])
    rel.project = scalar_projector(2,record['axes'])
    rel.pf_degree += 2
    assert rel.pf_degree > 13
    return rel


def rebuild_row(task):
    record,forms,sites = task
    rel = base_relation(record)
    assert [str(rel.cartan(w)) for w in sites] == record['sites']
    assert [str(sum(c*rel.coefficient(w) for w,c in form)) for form in forms] == record['scores']
    rel.clear()
    return True


def affected_sites(old, chosen, flavour):
    lookup = {w:i for i,w in enumerate(chosen)}
    result = []
    for w in old:
        entries = Counter()
        for i,t in enumerate(w):
            if t[:2] != (0,0) or not t[2+flavour]:
                continue
            seq = w[:i]+((1,1)+t[2:],)+w[i+1:]
            key = lambda v:(sum(v),v)
            odd = [key(v) for v in seq if parity(v)]
            if len(set(odd)) < len(odd):
                continue
            coordinate = tuple(sorted(seq,key=key))
            if coordinate in lookup:
                sign = (-1)**sum(a>b for j,a in enumerate(odd) for b in odd[j+1:])
                entries[lookup[coordinate]] += 2*t[2+flavour]*sign
        if entries:
            result.append((w,entries))
    return result


def run(workers,rebuild):
    old_rows,values,row = test(workers,rebuild)
    data = json.loads((HERE/'mixed.json').read_text())
    stored = data['rows']
    specs = list(mixed_specs())
    assert len(stored) == len(specs) == 76
    assert all((r['spec'],r['axes']) == pair for r,pair in zip(stored,specs,strict=True))
    small = fmpq_mat([r['sites'] for r in stored])
    coeffs = fmpq_mat([['19/64','1/4','1/4'],['5/6','1/3','1/3'],['3/4','1/12','1/12'],
                        ['5/96','-1/6','-1/6'],['1/4','-1/24','-1/24']])
    assert small*coeffs == fmpq_mat([r['scores'] for r in stored])
    pivots = [0,2,3,38,39]
    assert str(fmpq_mat([stored[i]['sites'] for i in pivots]).det()) == '19237055439568896000'
    forms = [ghost_form(row,x) for x in PHI]
    assert list(map(len,forms)) == [364,133,133]
    old = words((4,4,4,2,2))
    sites = [old[j] for j in data['plane_indices']]
    checks = range(76) if rebuild else pivots
    with ProcessPoolExecutor(max_workers=workers) as pool:
        assert all(pool.map(rebuild_row,[(stored[i],forms,sites) for i in checks]))
    new = words((5,5,4,2,2))
    fit = data['joint']
    cols = fit['cols']
    rows = [[r[j] for j in cols] for r in old_rows]
    ops = [affected_sites(old,[new[j] for j in cols],f) for f in range(3)]
    for i in pivots:
        record = stored[i]
        rel = base_relation(record)
        for f,terms in enumerate(ops):
            v = [Fraction() for _ in cols]
            for w,coordinates in terms:
                n = rel.commuting_coefficient(w)
                for j,m in coordinates.items():
                    v[j] += n*m
            rows.append(list(map(str,v)))
            values.append(record['scores'][f])
        rel.clear()
    mat = fmpq_mat(rows)
    assert mat.rank() == 23
    assert mat*fmpq_mat([[v] for v in fit['weights']]) == fmpq_mat([[v] for v in values])
    print('mixed base polarizations 76, re-evaluated',len(checks),'joint rational rank 23')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--rebuild',action='store_true')
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error('Run with assertions and a positive worker count')
    run(args.workers,args.rebuild)
