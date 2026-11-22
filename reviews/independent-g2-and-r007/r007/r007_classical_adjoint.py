"""Adjoint of the classical Taylor map in the free radial trace metric."""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path
from random import Random

import r007_free_radial as radial
import r007_primary_replay as jets


@lru_cache(None)
def vertex(u,v):
    p,q,h=u
    r,s,k=v
    if set(h).intersection(k):
        return None,0
    joined=p+r,q+s,tuple(sorted(h+k))
    coefficient=sum(c for c,w in jets.q_letter(joined) if w==(u,v))
    return joined, coefficient*radial.inner_letter(u)*radial.inner_letter(v)/radial.inner_letter(joined)


def glue(traces, marker):
    loc=[(i,j) for i,tr in enumerate(traces) for j,t in enumerate(tr) if t==marker]
    (a,p),(b,q)=loc
    if a==b:
        old=traces[a]
        return [tr for i,tr in enumerate(traces) if i!=a]+[old[p+1:q],old[q+1:]+old[:p]]
    left,right=traces[a],traces[b]
    return [tr for i,tr in enumerate(traces) if i not in (a,b)]+[left[p+1:]+left[:p]+right[q+1:]+right[:q]]


@lru_cache(None)
def sewn(lengths,i,j):
    # Slot 0 is the newly created letter. The two negative markers are
    # ordinary colour generators. Conjugation reverses their order.
    graph,offset = [[0,-2,-1]],0
    for length in lengths:
        graph.append([-1 if n==i else -2 if n==j else n+1 for n in range(offset,offset+length)])
        offset+=length
    graph=glue(glue(graph,-1),-2)
    if any(len(tr)==1 for tr in graph):
        return (),0
    return tuple(tuple(tr) for tr in graph if tr),sum(not tr for tr in graph)


def adjoint(state,n):
    result=Counter()
    for word,c in state.items():
        letters=tuple(x for tr in word for x in tr)
        lengths=tuple(map(len,word))
        parities=list(map(jets.odd,letters))
        counts=[sum(parities[:i]) for i in range(len(letters))]
        for i,u in enumerate(letters):
            for j,v in enumerate(letters):
                if i==j:
                    continue
                joined,coeff=vertex(u,v)
                if not coeff:
                    continue
                graph,loops=sewn(lengths,i,j)
                if not graph:
                    continue
                signed=parities[i]*counts[i]+parities[j]*(counts[j]-(i<j)*parities[i])
                indexed=(joined,)+letters
                order=[idx for tr in graph for idx in tr if jets.odd(indexed[idx])]
                signed+=sum(a>b for k,a in enumerate(order) for b in order[k+1:])
                new=tuple(tuple(indexed[idx] for idx in tr) for tr in graph)
                jets.add_term(result,new,(-1)**signed*n**loops*coeff*c)
    return jets.clean(result)


def trace(matrices):
    acc=matrices[0]
    for b in matrices[1:]:
        acc=[[sum(u*v for u,v in zip(row,col)) for col in zip(*b)] for row in acc]
    return sum(row[i] for i,row in enumerate(acc))


def matrix(rng,n):
    mat=[[rng.randrange(-2,3) for _ in range(n)] for _ in range(n)]
    mat[-1][-1]=-sum(mat[i][i] for i in range(n-1))
    return mat


def evaluate(normal,profile,n,seed):
    rng=Random(seed)
    even={x:matrix(rng,n) for x in sorted(set(profile)) if not jets.odd(x)}
    odd=sorted(x for x in profile if jets.odd(x))
    matrices=[matrix(rng,n) for x in odd]
    blocks={x:tuple(i for i,y in enumerate(odd) if x==y) for x in sorted(set(odd))}
    val=Fraction()
    for w,c in normal:
        for assign in product(*(permutations(l) for l in blocks.values())):
            byletter={x:iter(a) for x,a in zip(blocks,assign)}
            perm=[]
            factor=1
            for tr in w:
                selected=[]
                for x in tr:
                    if jets.odd(x):
                        p=next(byletter[x])
                        perm.append(p)
                        selected.append(matrices[p])
                    else:
                        selected.append(even[x])
                factor*=trace(selected)
            val+=c*factor*(-1)**sum(a>b for i,a in enumerate(perm) for b in perm[i+1:])
    return val


def run():
    state=radial.target()
    ans=adjoint(state,3)
    groups=defaultdict(list)
    for w,c in ans.items():
        profile=tuple(sorted(l for t in w for l in t))
        groups[profile].append((w,c))
    samples=[]
    for profile in sorted(groups,key=lambda g:(sum(map(jets.odd,g)),len(groups[g]),g)):
        val=evaluate(groups[profile],profile,3,596023)
        if val:
            samples.append(dict(profile=profile,terms=len(groups[profile]),seed=596023,value=str(val)))
            if len(samples)==3:
                break
    return dict(target_terms=len(state),adjoint_trace_terms=len(ans),adjoint_profiles=len(groups),
                nonzero_finite_SU3_evaluations=samples)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    if not __debug__:
        raise RuntimeError('run with assertions enabled')
    record=run()
    a.out.write_text(json.dumps(record,indent=2)+'\n')
    print(record,flush=True)
