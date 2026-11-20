"""Additional exceptional-colour boundaries on the F4 second- and third-order source spaces."""
import argparse
import math
import random
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from f4_colour import basis, frame, jordan_colours
from f4_formal_cycles import coefficient_matrix, probe
from f4_joined_boundaries import evaluate, graph, native
from f4_singlet_characters import singlets
from fortuity.full.letters import Q_letter, is_ferm
from mixed33_exact import read
from mixed_middle import pin


def inversions(seq):
    return sum(a>b for i,a in enumerate(seq) for b in seq[i+1:])


@lru_cache(None)
def short_arm(word):
    if len(word)!=4:
        return 1,word
    odd=tuple(x for x in word if is_ferm(x))
    if len(set(odd))!=len(odd):
        return 0,word
    return (-1)**inversions(odd),tuple(sorted(word))


@lru_cache(maxsize=65536)
def canonical(pair):
    (s,a),(t,b)=map(short_arm,pair)
    if not s*t:
        return 0,pair
    sign,g=graph((a,b),"jordan")
    return s*t*sign,g


def check_short_arms():
    matrices,_=basis()
    jordan=jordan_colours()
    assert not np.einsum("iab,jba->ij",jordan,matrices).any()
    pairs=matrices[:,None]@matrices[None,:]
    assert not np.einsum("iab,jkba->ijk",jordan,pairs).any()
    # Interchanging two Lie generators in a cubic inserts one bracket, leaving
    # a quadratic word. After symmetrization, transposition kills the cubic.
    return [26*52,26*52**2]


@lru_cache(None)
def spin_flavours(t,f):
    spins=[]
    for bits in product((0,1),repeat=t):
        spins.append(((-1)**sum(bits),tuple(d for x in bits for d in (x,1-x))))
    flavours=[]
    for copies in product(tuple(permutations(range(3))),repeat=f):
        flavours.append(((-1)**sum(map(inversions,copies)),tuple(i for p in copies for i in p)))
    return [(s*v,a,b) for s,a in spins for v,b in flavours]


@lru_cache(maxsize=100000)
def slot_letter(spec,spin,flavour,spin_slots):
    z=[spin[i] for i in spec if i<spin_slots]
    th=[flavour[i-spin_slots] for i in spec if i>=spin_slots]
    if len(set(th))!=len(th):
        return 0,None
    a,b=z.count(0),z.count(1)
    c=math.factorial(a)*math.factorial(b)*(-1)**inversions(th)
    return c,(a,b,tuple(sorted(th)))


def source(t,f,deficit,seed):
    rng=random.Random(seed)
    slots=list(range(2*t+3*f))
    rng.shuffle(slots)
    if deficit==1:
        cuts=(2,)
    elif deficit==2:
        cuts=(2,2) if seed%2 else (3,)
    else:
        raise ValueError(deficit)
    pattern=[]
    offset=0
    for length in cuts:
        pattern.append(tuple(slots[offset:offset+length]))
        offset+=length
    pattern.extend((x,) for x in slots[offset:])
    rng.shuffle(pattern)
    n=len(pattern)
    cut=rng.randrange(4,n-3)
    result=Counter()
    for sign,spin,flavours in spin_flavours(t,f):
        w=[]
        for cell in pattern:
            c,l=slot_letter(cell,spin,flavours,2*t)
            sign*=c
            if not sign:
                break
            w.append(l)
        if sign:
            sg,g=canonical((tuple(w[:cut]),tuple(w[cut:])))
            if sg:
                result[g]+=sg*sign
    return {g:c for g,c in result.items() if c}


def variation(expr,axis):
    result=Counter()
    for g,c in expr.items():
        phase=1
        for j,w in enumerate(g):
            for k,l in enumerate(w):
                a,b,theta=l
                if axis==0:
                    replacements=[(phase*n,(x,y)) for n,x,y in Q_letter(l)]
                elif axis==3:
                    replacements=[(a+1,((a+1,b-1,theta),))] if b else []
                else:
                    flavour=tuple(sorted(axis-1 if x==axis else x for x in theta))
                    replacements=[(1,((a,b,flavour),))] if axis in theta and axis-1 not in theta else []
                for n,letters in replacements:
                    u=w[:k]+letters+w[k+1:]
                    s,v=canonical(g[:j]+(u,)+g[j+1:])
                    if s:
                        result[v]+=s*c*n
                if is_ferm(l):
                    phase=-phase
    return {g:c for g,c in result.items() if c}


def run(t,kind,plan,previous,work,count,workers):
    work.mkdir(exist_ok=True,parents=True)
    z,old=read(plan),read(previous)
    record=next(e for e in old["records"] if e["weight"][:2]==[t,t] and e["S_factors"]==1)
    assert pin(plan)==record[kind]["plan_sha256"]
    p=old["prime"]
    weight=tuple(z["mu"])
    assert weight==(t,t)+(4-t,)*3
    deficiency=sum(weight)-z["length"]
    assert deficiency==(2 if kind=="incoming" else 1)
    words=list(map(m.decode,z["target_dictionary"]))
    dims=singlets("F4",weight,z["length"]+1)
    assert set(dims)=={m.nferm(w) for w in words}
    saved=record[kind]["rank"]
    rows=[z["images"][i] for i in saved["rows"]]
    base,frames=probe(rows,words,dims,p,m.library(work,27,p),workers,1507032+z["length"]*229)
    assert int(nmod_mat(base[:,saved["columns"]].tolist(),p).det())==saved["determinant"]
    letters={d:{l for w in words if m.nferm(w)==d for tr in w for l in tr} for d in dims}
    lib=native(work,p,"jordan")
    fframes=[frame(letters[x["degree"]],x["points"],x["k"],x["degree"],p,x["seed"]) for x in frames]
    seeds=[];matrices=[base];graph_counts=[]
    checked=check_short_arms()
    start=329404+t*10000+deficiency*1100
    for first in range(start,start+count,80):
        recipes=[]
        for seed in range(first,min(first+80,start+count)):
            expr=source(t,4-t,deficiency,seed)
            if not expr:
                continue
            assert all(not variation(expr,axis) for axis in (1,2,3))
            varied=variation(expr,0)
            if varied:
                assert not variation(varied,0)
                seeds.append(seed)
                recipes.append(varied)
        if not recipes:
            continue
        blocks=[];total=0
        for fr,x in zip(fframes,frames):
            d=x["degree"]
            gs=sorted({g for r in recipes for g in r if m.nferm(g)==d})
            if not gs:
                blocks.append(np.zeros((len(recipes),fr.ncols),np.int64))
                continue
            assert all(l in letters[d] for g in gs for tr in g for l in tr)
            ids={g:i for i,g in enumerate(gs)}
            co=coefficient_matrix([[(ids[g],v) for g,v in r.items() if g in ids] for r in recipes],len(gs),p)
            blocks.append(np.asarray(co@evaluate(gs,fr,work,workers,lib,"jordan"))%p)
            total+=len(gs)
        graph_counts.append(total)
        matrices.append(np.column_stack(blocks))
        minor=m.minor(np.vstack(matrices),p)
        result=dict(prime=p,weight=weight,kind=kind,deficit=deficiency,plan_sha256=pin(plan),
                    baseline=saved["rank"],short_arm_checks=checked,attempted=min(first+80,start+count)-start,
                    first_seed=start,seeds=seeds,frame_specs=frames,graphs_by_batch=graph_counts,minor=minor)
        m.write(work/(kind+".json"),result)
        print(t,kind,"parents",len(seeds),"rank",minor["rank"],"previous",saved["rank"],flush=True)
        if minor["rank"]==sum(dims.values()):
            break
    if len(matrices)==1:
        m.write(work/(kind+".json"),dict(prime=p,kind=kind,weight=weight,plan_sha256=pin(plan),
                baseline=saved["rank"],short_arm_checks=checked,attempted=count,first_seed=start,
                seeds=[],minor=saved))


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    for key in ("plan","previous","work"):
        a.add_argument("--"+key,type=Path,required=True)
    a.add_argument("--t",choices=(1,2),type=int,required=True)
    a.add_argument("--kind",choices=("incoming","outgoing"),required=True)
    a.add_argument("--count",type=int,default=500)
    a.add_argument("--workers",type=int,default=12)
    args=a.parse_args()
    assert __debug__
    run(**vars(args))
