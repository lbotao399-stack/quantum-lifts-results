"""Convert protected highest-weight trace products to native classical-jet coordinates."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import factorial
import gzip
import json
from pathlib import Path
import sys


def native(letter):
    plus, minus, flavours = letter
    degree = len(flavours)
    if degree == 0:
        typ, flavour, sign = 3, 0, -1
    elif degree == 1:
        typ, flavour, sign = 0, flavours[0], 1
    elif degree == 2:
        missing, = set(range(3)) - set(flavours)
        typ, flavour, sign = 1, missing, (-1)**missing
    else:
        typ, flavour, sign = 2, 0, 1
    return (typ,flavour,plus,minus), Fraction(sign,factorial(plus)*factorial(minus))


def check_dictionary(max_degree):
    from fortuity.full.letters import Q_letter
    from math import comb
    count = 0
    for plus in range(max_degree+1):
        for minus in range(max_degree+1-plus):
            for degree in range(4):
                if plus+minus+degree == 0 or plus+minus+degree>max_degree:
                    continue
                for flavours in combinations(range(3),degree):
                    (typ,fl,p,q), coefficient = native((plus,minus,flavours))
                    components = Counter()
                    for c, left, right in Q_letter((p,q,flavours)):
                        a,x = native(left)
                        b,y = native(right)
                        components[a,b] += c*x*y
                    expected = Counter()
                    def bracket(t,f,u,g,scale):
                        for i in range(p+1):
                            for j in range(q+1):
                                left,right=(t,f,i,j),(u,g,p-i,q-j)
                                if (t==3 and i+j==0) or (u==3 and p-i+q-j==0):
                                    continue
                                c=scale*coefficient*comb(p,i)*comb(q,j)
                                expected[left,right]+=c
                                expected[right,left]+=c if t in (1,3) and u in (1,3) else -c
                    if typ==3:
                        bracket(3,0,3,0,Fraction(-1,2))
                    else:
                        bracket(3,0,typ,fl,-1)
                        if typ==1:
                            bracket(0,(fl+1)%3,0,(fl+2)%3,1)
                        if typ==2:
                            for i in range(3):
                                bracket(0,i,1,i,-1)
                    assert components == expected,(p,q,flavours,components,expected)
                    count+=1
    return count


def run(args):
    sys.path.insert(0,str(args.repo.resolve()/'algorithms/dual_singlets_20261118'))
    import mixed_incoming as m
    from sp3_omega_product_replay import trace_factor
    from sp3_omega_products import raising
    packet=json.loads(gzip.decompress(args.recipes.read_bytes()))
    combinations_=[[([tuple(g) for g in spec],c) for spec,c in row] for row in packet['recipes']]
    assert check_dictionary(sum(packet['weight'])-packet['length']+1)>0
    for row in combinations_:
        raised=Counter()
        for spec,c in row:
            for w,d in raising(tuple(spec)).items():
                raised[w]+=c*d
        assert not any(raised.values())
    all_factors={g for row in combinations_ for spec,c in row for g in spec}
    assert all(not m.q(trace_factor(g)) for g in all_factors)
    from functools import lru_cache
    @lru_cache(None)
    def on_trace(trace):
        out,scale=[],Fraction(1)
        for jet in trace:
            l,c=native(jet)
            out.append(l)
            scale*=c
        return tuple(out),scale
    terminal={}
    traces={}
    row_terms=[]
    @lru_cache(maxsize=5000)
    def word_product(spec):
        terms = {(): Fraction(1)}
        for factor in spec:
            out=Counter()
            for w,c in terms.items():
                for t,d in trace_factor(factor).items():
                    s,k=m.so_canon(w+t)
                    if s: out[k]+=s*c*d
            terms={k:c for k,c in out.items() if c}
        return terms
    for row in combinations_:
        terms=Counter()
        for spec,c in row:
            for word,d in word_product(tuple(spec)).items():
                terms[word]+=c*d
        sparse=[]
        for w,c in terms.items():
            if not c: continue
            ids=[]
            for t in w:
                l,s=on_trace(t)
                ids.append(traces.setdefault(l,len(traces)))
                c*=s
            sparse.append((terminal.setdefault(tuple(ids),len(terminal)),c))
        row_terms.append(sorted(sparse))
    with args.out.open('w') as out:
        print(len(traces),len(terminal),len(row_terms),file=out)
        for trace in traces:
            print(len(trace),*(v for jet in trace for v in jet),file=out)
        for w in terminal:
            print(len(w),*w,file=out)
        for row in row_terms:
            print(len(row),*(f'{i} {c.numerator} {c.denominator}' for i,c in row),file=out)
    print('traces, words, rows, coefficients',len(traces),len(terminal),len(row_terms),sum(map(len,row_terms)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--recipes',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        raise RuntimeError('enable Python assertions')
    run(a)
