"""Check the odd substitutions on every allowed local letter and ordered input pair."""
import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
import sys


def run(a):
    sys.path.insert(0, str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from rank_four_lowered_replay import Q1d_multi
    from q1d_exact import pair_terms
    from mixed_incoming import so_canon
    from so8_level24_descendant import check_letter, substitute, substitution
    from mixed_ward import transform
    jets = {n: [(i,j,t) for i in range(n+1) for j in range(n+1-i)
                for t in combinations(range(3), n-i-j)] for n in (1,2)}
    assert [len(jets[n]) for n in (1,2)] == [5,12]
    steps = list(product(range(2), range(3)))
    alphabet = jets[1]+jets[2]
    for step in steps:
        for l in alphabet:
            check_letter(l,*step,True)
        for s in steps:
            for l in alphabet:
                row = Counter()
                for x,y in ((s,step),(step,s)):
                    v,m = substitution(l,*x,True)
                    if v:
                        n,k = substitution(m,*y,True)
                        if n:
                            row[k] += v*n
                assert not any(row.values())
    ordered = list(product(jets[1], jets[1]))+list(product(jets[1],jets[2]))+list(product(jets[2],jets[1]))
    active = [pair for pair in ordered if pair_terms(*pair)]
    for left,right in ordered:
        for step in steps:
            n,u = substitution(left,*step,True)
            m,v = substitution(right,*step,True)
            if n:
                assert not pair_terms(u,right)
            if m:
                assert not pair_terms(left,v)
            for c,l,r in pair_terms(left,right):
                assert not l[2] and not r[2]
                assert not substitution(l,*step,True)[0] and not substitution(r,*step,True)[0]
    # Check disjoint odd actions with genuine joins, splits and orthogonal reversals.
    testing = set()
    xx = jets[1]
    for k,(l,r) in enumerate(active):
        for i,u in enumerate(xx):
            for j,v in enumerate(xx):
                other = xx[(i+2*j+k)%len(xx)]
                for traces in (((l,u,r,v,other),), ((u,l,v),(r,other)), ((u,r,l),(other,v))):
                    sign,w = so_canon(traces)
                    if sign:
                        testing.add(w)
    def q1(expr):
        result = Counter()
        for w,c in expr.items():
            for raw,v in Q1d_multi(w,13,'SO_odd').items():
                s,k = so_canon(raw)
                if s:
                    result[k] += s*c*v
        return {w:v for w,v in result.items() if v}
    for w in testing:
        direct = q1({w:1})
        for step in steps:
            row = Counter(substitute(direct,*step,True))
            row.update(q1(substitute({w:1},*step,True)))
            assert not any(row.values()), (w,step,row)
    g = (1,1,0,0,0)
    lowered = transform(transform({(g,)*7:1},0,0,False),1,0,False)
    basis = [tuple(map(tuple,row)) for row in json.loads(a.case.read_text())['target']]
    coordinates = [lowered.get(w,0) for w in basis]
    assert len(lowered)==2 and coordinates[0] and 6*coordinates[0]==coordinates[1] and not coordinates[2]
    result = dict(first_jets=len(jets[1]), second_jets=len(jets[2]),
                  classical_anticommutators=len(alphabet)*len(steps),
                  substitution_anticommutators=len(alphabet)*len(steps)**2,
                  ordered_pairs=len(ordered), active_ordered_pairs=len(active),
                  ordered_pair_substitutions=len(ordered)*len(steps),
                  native_multitraces=len(testing), native_trace_equations=len(testing)*len(steps),
                  twice_lowered_S7=coordinates)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(result,flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('repo','case','out'):
        p.add_argument('--'+k,type=Path,required=True)
    assert __debug__
    run(p.parse_args())
