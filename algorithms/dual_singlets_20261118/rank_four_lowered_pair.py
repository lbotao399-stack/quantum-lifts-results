"""First quantum action after the six classical spin-to-flavour substitutions on the SO(9) source."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from itertools import accumulate
from pathlib import Path

import mixed_incoming as m
from fortuity.full.letters import Q_letter, is_ferm
from mixed33_exact import read
from rank_four_s5 import Q1d_multi
from q1d_exact import pair_terms
from rank_four_sixth import unpack
from replay_orthogonal import exchanges


SIMPLE = ((1,0,()),(0,1,()),(0,0,(0,)),(0,0,(1,)),(0,0,(2,)))
STEPS = tuple((spin,flavour) for spin in range(2) for flavour in range(3))


@lru_cache(None)
def replace(letter, spin, flavour):
    a,b,theta = letter
    if not letter[spin] or flavour in theta:
        return None
    sign = (-1)**sum(i > flavour for i in theta)
    return sign,(a-(spin==0),b-(spin==1),tuple(sorted(theta+(flavour,))))


def lower(letters, spin, flavour):
    sign = 1
    for j,letter in enumerate(letters):
        v = replace(letter,spin,flavour)
        if v:
            n,new = v
            yield sign*n,letters[:j]+(new,)+letters[j+1:]
        if is_ferm(letter):
            sign = -sign


def check_classical(alphabet):
    zero = 0
    for letter in alphabet:
        for spin,flavour in STEPS:
            residual = Counter()
            for c,(new,) in lower((letter,),spin,flavour):
                for n,a,b in Q_letter(new):
                    residual[a,b] += c*n
            for n,a,b in Q_letter(letter):
                for c,pair in lower((a,b),spin,flavour):
                    residual[pair] += c*n
            assert not any(residual.values())
            zero += 1
    return zero


@lru_cache(None)
def routings(shape,i,j,least):
    offsets = [0,*accumulate(shape)]
    traces = [tuple(range(lo,hi)) for lo,hi in zip(offsets,offsets[1:])]
    p,q = [next(k for k,tr in enumerate(traces) if idx in tr) for idx in (i,j)]
    rest = tuple(tr for k,tr in enumerate(traces) if k not in (p,q))
    out = []
    for x in (0,1):
        for y in (0,1):
            sign = (-1)**(x+y)
            u,v = i-offsets[p]+x,j-offsets[q]+y
            if p != q:
                a,b = traces[p][:u],traces[p][u:]
                d,e = traces[q][:v],traces[q][v:]
                out.append((sign,rest+(a+e+d+b,)))
                out.append((-sign*(-1)**len(traces[q]),rest+(a+d[::-1]+e[::-1]+b,)))
            else:
                lo,hi = sorted((u,v))
                a,b,c = traces[p][:lo],traces[p][lo:hi],traces[p][hi:]
                out.append((sign*9**(int(not b)+int(not (a+c))),
                            rest+tuple(t for t in (b,a+c) if t)))
                out.append((-sign*(-1)**len(b),rest+(a+b[::-1]+c,)))
    return [(c,route) for c,route in out if len(route)>=least and all(len(tr)%2==0 for tr in route)]


def diagrams(word, least):
    shape = tuple(map(len,word))
    letters = sum(word,())
    double, = [j for j,l in enumerate(letters) if sum(l[:2])+len(l[2])==2]
    result = Counter()
    for other in range(len(letters)):
        if other == double:
            continue
        i,j = sorted((double,other))
        typ = [-2 if k==i else -1 if k==j else SIMPLE.index(l) for k,l in enumerate(letters)]
        for c,route in routings(shape,i,j,least):
            blocks = [tuple(sorted(tr,key=typ.__getitem__)) for tr in route]
            def labels(tr):
                return tuple(typ[k] for k in tr)
            desc = [labels(tr) for tr in blocks]
            if any(d.count(0)>1 or d.count(1)>1 for d in desc):
                continue
            blocks.sort(key=lambda tr:(len(tr),labels(tr)))
            desc = tuple(labels(tr) for tr in blocks)
            if any(u==v and (u.count(0)+u.count(1))%2 for u,v in zip(desc,desc[1:])):
                continue
            order = (i,j)+tuple(k for tr in blocks for k in tr if k not in (i,j))
            phase = (-1)**exchanges([k for k in order if is_ferm(letters[k])])
            result[(letters[i],letters[j],desc)] += phase*c
    return {key:c for key,c in result.items() if c}


def input_and_blocks(key):
    a,b,desc = key
    letters = [a,b]
    blocks = []
    for tr in desc:
        ids = []
        for t in tr:
            if t < 0:
                ids.append(t+2)
            else:
                ids.append(len(letters))
                letters.append(SIMPLE[t])
        blocks.append(tuple(ids))
    return tuple(letters),tuple(blocks)


def reduced(species,blocks):
    pieces = []
    for tr in blocks:
        seq = sorted(tr,key=species.__getitem__)
        counts = tuple(sum(species[k]==s for k in seq) for s in range(5))
        if counts[0]>1 or counts[1]>1:
            return 0,()
        pieces.append((counts,seq))
    pieces.sort(key=lambda x:(sum(x[0]),x[0]))
    if any(a==b and (a[0]+a[1])%2 for (a,_),(b,_) in zip(pieces,pieces[1:])):
        return 0,()
    phase = (-1)**exchanges([i for _,seq in pieces for i in seq if species[i]<2])
    return phase,tuple(w for w,_ in pieces)


def projected_pair(letters,blocks,lookup):
    value = [0]*len(next(iter(lookup.values())))
    for c,a,b in pair_terms(*letters[:2]):
        full = (a,b)+letters[2:]
        if any(l not in SIMPLE for l in full):
            continue
        phase,word = reduced(tuple(SIMPLE.index(l) for l in full),blocks)
        if phase and word in lookup:
            pair = lookup[word]
            for j,n in enumerate(pair):
                value[j] += phase*c*n
    return value


def evaluate(task):
    key, lookup, steps = task
    letters,blocks = input_and_blocks(key)
    states = {letters:1}
    for spin,flavour in steps:
        new = Counter()
        for word,v in states.items():
            for n,w in lower(word,spin,flavour):
                new[w] += v*n
        states = {w:v for w,v in new.items() if v}
    total = [0]*len(next(iter(lookup.values())))
    for w,n in states.items():
        if len(w[0][2])+len(w[1][2])!=3:
            continue
        values = projected_pair(w,blocks,lookup)
        for j,c in enumerate(values):
            total[j] += n*c
    return total


def direct(word,lookup):
    values = [0]*len(next(iter(lookup.values())))
    for term,c in Q1d_multi(word,9,group="SO_odd").items():
        if any(len(tr)%2 for tr in term):
            continue
        species = tuple(SIMPLE.index(l) for tr in term for l in tr)
        ends = [0,*accumulate(map(len,term))]
        phase,w = reduced(species,[range(lo,hi) for lo,hi in zip(ends,ends[1:])])
        if phase and w in lookup:
            for j,n in enumerate(lookup[w]):
                values[j] += c*phase*n
    return values


def controls(lookup,least):
    import random
    rng = random.Random(67119)
    tested = nonzero = 0
    for trial in range(240):
        shape = ((2,2,2,2,2,2),(2,2,2,2,4),(2,2,4,4),(2,3,3,4))[trial%4]
        f = trial%3
        letters = [SIMPLE[0]]*2+[SIMPLE[1]]*2+[(0,0,tuple(i for i in range(3) if i!=f))]
        letters.extend(SIMPLE[2+i] for i in range(3) for _ in range(2+int(i==f)))
        rng.shuffle(letters)
        off = [0,*accumulate(shape)]
        word = tuple(tuple(letters[a:b]) for a,b in zip(off,off[1:]))
        expected = direct(word,lookup)
        value = [0]*len(next(iter(lookup.values())))
        for key,n in diagrams(word,least).items():
            vv = evaluate((key,lookup,()))
            for j,c in enumerate(vv):
                value[j] += n*c
        assert value == expected,(word,value,expected)
        tested += 1
        nonzero += int(any(value))
    assert nonzero
    return tested,nonzero


MULTIPLIERS = {
    "gauginos": ((SIMPLE[0],SIMPLE[1]),STEPS),
    "mixed": ((SIMPLE[0],SIMPLE[2]),tuple(t for t in STEPS if t != (1,0))),
    "scalars": ((SIMPLE[2],SIMPLE[2]),tuple(t for t in STEPS if t[1] != 0)),
}


def run(source, classes, out, workers, prepare=False, channel="gauginos", cartan=None):
    cert = read(classes)
    basis = [tuple(map(tuple,w)) for w in cert["formal_weight_basis"]]
    columns = cert["weight_functionals"]
    if cartan is not None:
        from flint import fmpz_mat
        from math import gcd
        from rank_four_terminal import quotient
        sp = read(cartan)
        basis_again,cov = quotient(tuple(cert["weight"]),10)
        assert basis_again == tuple(basis) or basis_again == basis
        relations = fmpz_mat([r["singlet"] for r in sp["relations"]])
        kernel,dim = relations.nullspace()
        assert dim==3
        for k in range(dim):
            v = fmpz_mat([[kernel[i,k]] for i in range(cov.ncols())])
            assert (relations*v).is_zero()
            coeff = list(map(int,(cov*v).entries()))
            den = gcd(*coeff)
            columns.append([[i,str(n//den)] for i,n in enumerate(coeff) if n])
    lookup = {}
    for j,items in enumerate(columns):
        for i,n in items:
            lookup.setdefault(basis[i],[0]*len(columns))[j] = Fraction(n)
    least = min(map(len,lookup))
    ctrl = controls(lookup,least)
    bilinear, steps = MULTIPLIERS[channel]
    src = Counter()
    for w,n in unpack(source).items():
        sign,word = m.so_canon((bilinear,)+w)
        if sign:
            src[word] += sign*n
    src = {w:n for w,n in src.items() if n}
    alphabet = {l for w in src for tr in w for l in tr}
    checked = check_classical(alphabet)
    combined = Counter()
    for w,c in src.items():
        if len(w)+1 >= least:
            for key,n in diagrams(w,least).items():
                combined[key] += n*c
    keys = [(key,n) for key,n in combined.items() if n]
    print("source",len(src),"classical checks",checked,"projected diagrams",len(keys),
          "direct controls",ctrl,flush=True)
    if prepare:
        return
    total = [0]*len(columns)
    with ProcessPoolExecutor(workers) as pool:
        for (_,n),v in zip(keys,pool.map(evaluate,((key,lookup,steps) for key,_ in keys),chunksize=2)):
            for j,c in enumerate(v):
                total[j] += n*c
    m.write(out,dict(group="SO(9)",level=30,source="spin-to-flavour actions on a bilinear times F",
                     bilinear=bilinear,steps=steps,input_terms=len(src),projected_diagrams=len(keys),
                     classical_anticommutators=checked,direct_controls=ctrl,
                     projected_first_action=list(map(str,total))))
    print("first action",total,flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("source","classes","out"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--prepare",action="store_true")
    p.add_argument("--channel",choices=MULTIPLIERS,default="gauginos")
    p.add_argument("--cartan",type=Path)
    assert __debug__
    run(**vars(p.parse_args()))
