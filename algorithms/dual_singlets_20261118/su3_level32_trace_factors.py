"""Test the level-32 first image on a subspace annihilating quadratic traces."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import permutations, product
import json
from pathlib import Path
from random import Random

from su3_level32_replay import DATA
from block_first import code, odd


U = (-1,1,0,-1,0,0,0,0,1)
V = (0,0,0,0,0,1,0,0,0)
W = (0,0,0,0,0,0,1,0,0)
IDENTITY = (1,0,0,0,1,0,0,0,1)


def multiply(a,b):
    return tuple(a[i]*b[j]+a[i+1]*b[3+j]+a[i+2]*b[6+j]
                 for i in (0,3,6) for j in (0,1,2))


def trace(a):
    return a[0]+a[4]+a[8]


def linear(matrices,coefficients):
    return tuple(sum(c*m[k] for c,m in zip(coefficients,matrices,strict=True)) for k in range(9))


def evaluator(t, basis, seed):
    rng = Random(seed)
    odd_counts = Counter(k for k in t if odd(k))
    matrices = {}
    exterior = {}
    start = 0
    for k in sorted(set(t)):
        indices = list(range(start,start+odd_counts[k])) if odd(k) else [-1]
        if odd(k):
            start += len(indices)
            exterior[k] = indices
        for i in indices:
            matrices[k,i] = linear(basis,[rng.randrange(-3,4) for _ in basis])

    @lru_cache(30000)
    def tr(word,labels):
        indices = iter(labels)
        matrix = IDENTITY
        for k in word:
            matrix = multiply(matrix,matrices[k,next(indices) if odd(k) else -1])
        return trace(matrix)

    assignments=list(product(*(tuple(permutations(v)) for v in exterior.values())))
    @lru_cache(40000)
    def evaluate(w):
        ev=0
        for choice in assignments:
            indices={k:iter(v) for k,v in zip(exterior,choice,strict=True)}
            exterior_order=[]
            term=1
            for word in w:
                sub=tuple(next(indices[k]) for k in word if odd(k))
                term *= tr(word,sub)
                if not term:
                    break
                exterior_order.extend(sub)
            if term:
                ev+=(-1)**sum(k>l for i,k in enumerate(exterior_order) for l in exterior_order[i+1:])*term
        return ev
    return evaluate,dict(exterior_degree=start,assignments=[[k,i,m] for (k,i),m in sorted(matrices.items())])


def find(data, basis, seed):
    gram=[[trace(multiply(a,b)) for b in basis] for a in basis]
    cartan=[trace(multiply(multiply(a,b),c))-trace(multiply(multiply(b,a),c))
            for a,b,c in product(basis,repeat=3)]
    assert all(trace(a)==0 for a in basis) and not any(map(any,gram))
    rows=defaultdict(list)
    counts=Counter()
    h=hashlib.sha256()
    with gzip.open(data/"Q1.jsonl.gz","rb") as stream:
        for line in stream:
            h.update(line)
            c,raw=json.loads(line)
            if any(len(tr)==2 for tr in raw):
                continue
            w=tuple(bytes(map(code,tr)) for tr in raw)
            t=bytes(sorted(b"".join(w)))
            odd_counts=Counter(i for i in t if odd(i))
            if max(odd_counts.values(),default=0)>len(basis):
                continue
            rows[t].append((w,Fraction(c)))
            counts[t]+=1
    assert h.hexdigest()==json.loads((data/"order1.json").read_bytes())["output_sha256"]
    checked=0
    for t in sorted(rows,key=lambda t:(counts[t],t)):
        checked+=1
        evaluate, point=evaluator(t,basis,seed)
        numerator=sum((c*evaluate(w) for w,c in rows[t]),Fraction())
        if numerator:
            return dict(source_sha256=h.hexdigest(),colour_subspace=basis,quadratic_gram=gram,
                        traces_of_commutators_times_basis=cartan,seed=seed,contents_examined=checked,
                        selected_content=list(t),retained_words=counts[t],Q1_source=str(numerator),**point)
    return dict(source_sha256=h.hexdigest(),colour_subspace=basis,contents_examined=checked,
                seed=seed,Q1_source="0")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--seed",type=int,default=158399)
    p.add_argument("--dimension",choices=(2,3),type=int,default=2)
    p.add_argument("--commutator-orthogonal",action="store_true")
    a=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    basis=(linear((U,V,W),(1,2,-1)),linear((U,V,W),(1,-1,2))) if a.dimension==2 else (U,V,W)
    if a.commutator_orthogonal:
        basis=((0,0,1,0,0,0,0,0,0),(0,1,0,0,0,1,0,0,0),(0,0,0,1,0,0,0,-1,0))
    result=find(a.data,basis,a.seed)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(a.dimension,result["contents_examined"],result["Q1_source"])
