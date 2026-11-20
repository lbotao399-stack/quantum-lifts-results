"""Replay the two independent SO(9) first actions by transforming whole canonical trace words."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import mixed_incoming as m
from fortuity.full.letters import is_ferm
from mixed33_exact import read
from rank_four_s5 import Q1d_multi
from rank_four_sixth import unpack
from orthogonal_relations import norm


def derivative(word, spin, flavour):
    result = Counter()
    sign = 1
    for i,tr in enumerate(word):
        for j,(a,b,theta) in enumerate(tr):
            if (a,b)[spin] and flavour not in theta:
                # The inserted derivative is ordered after the existing flavour derivatives.
                c = (-1)**sum(k>flavour for k in theta)
                new = (a-int(spin==0), b-int(spin==1), tuple(sorted((*theta,flavour))))
                normalized,term = m.so_canon(word[:i]+(tr[:j]+(new,)+tr[j+1:],)+word[i+1:])
                if normalized:
                    result[term] += sign*c*normalized
            if is_ferm((a,b,theta)):
                sign = -sign
    return {w:n for w,n in result.items() if n}


def terminal(word):
    blocks = []
    sign = 1
    for tr in word:
        if len(tr)%2:
            return 0,()
        bosons = [0]*3
        odd = []
        for a,b,theta in tr:
            if theta:
                assert a==b==0 and len(theta)==1
                bosons[theta[0]] += 1
            else:
                assert a+b==1
                odd.append(0 if a else 1)
        if len(set(odd))!=len(odd):
            return 0,()
        if odd==[1,0]:
            sign = -sign
        blocks.append((odd.count(0),odd.count(1),*bosons))
    c,w = norm(tuple(blocks))
    return c*sign,w


@lru_cache(16384)
def direct_action(word):
    values = [0,0]
    for x,c in Q1d_multi(word,9,group="SO_odd").items():
        if len(x)<5:
            continue
        sign,monomial = terminal(x)
        if sign and monomial in LOOKUP:
            for j,n in enumerate(LOOKUP[monomial]):
                values[j] += c*sign*n
    return tuple(values)


def init_worker(lookup):
    global LOOKUP
    LOOKUP = lookup


def evaluate(task):
    word,steps = task
    expr = {word:1}
    for spin,flavour in steps:
        terms = Counter()
        for old,n in expr.items():
            for w,c in derivative(old,spin,flavour).items():
                terms[w] += n*c
        expr = {w:n for w,n in terms.items() if n}
    values = [0,0]
    for w,n in expr.items():
        if not any(len(l[2])==2 for tr in w for l in tr):
            continue
        for j,c in enumerate(direct_action(w)):
            values[j] += n*c
    return values


def run(source, classes, output, workers):
    cert = read(classes)
    basis = [tuple(map(tuple,w)) for w in cert["formal_weight_basis"]]
    lookup = {}
    for j,column in enumerate(cert["weight_functionals"]):
        for i,n in column:
            lookup.setdefault(basis[i],[0,0])[j] = Fraction(n)
    assert min(map(len,lookup))==5
    src = unpack(source)
    result = {}
    with ProcessPoolExecutor(workers,initializer=init_worker,initargs=(lookup,)) as pool:
        for name,bilinear,flavours in (
                ("gauginos",((1,0,()),(0,1,())),(0,1,2)),
                ("scalars",((0,0,(0,)),(0,0,(0,))),(1,2))):
            steps = tuple((spin,flavour) for spin in (0,1) for flavour in flavours)
            selected = []
            for old,n in src.items():
                sign,w = m.so_canon((bilinear,)+old)
                if not sign or len(w)<4 or any(len(tr)%2 for tr in w):
                    continue
                selected.append((w,n*sign))
            cols = [0,0]
            stream = pool.map(evaluate,((w,steps) for w,n in selected),chunksize=4)
            for (_,n),values in zip(selected,stream):
                for j,c in enumerate(values):
                    cols[j] += n*c
            result[name] = dict(source_words=len(selected),projected_first_action=list(map(str,cols)))
            print(name,len(selected),result[name]["projected_first_action"],flush=True)
    a,b = (list(map(Fraction,result[k]["projected_first_action"])) for k in ("gauginos","scalars"))
    m.write(output,dict(channels=result,determinant=str(a[0]*b[1]-a[1]*b[0])))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("source","classes","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    assert __debug__
    run(**vars(p.parse_args()))
