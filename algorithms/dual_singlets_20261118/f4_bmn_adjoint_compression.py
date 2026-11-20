"""Classical BMN boundaries from cyclic products of Albert words projected to the adjoint."""
import argparse
import hashlib
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

import mixed_incoming as m
import scalar_neighbors as scalar
from f4_bmn_double_jordan import inherited
from f4_bmn_primary import LETTERS, SIGNS, WEIGHT, bmn_convention
from f4_colour import basis, frame
from f4_formal_cycles import coefficient_matrix
from f4_joined_boundaries import native
from mixed33_exact import read


@lru_cache(160000)
def canonical(blocks):
    degrees = [sum(3<=v<6 for v in word) for word in blocks]
    total = sum(degrees)
    transposition = (-1)**(sum(map(len,blocks))+total*(total-1)//2)
    choices = []
    for words,counts,sign in ((blocks,degrees,1),(tuple(w[::-1] for w in blocks[::-1]),degrees[::-1],transposition)):
        before=0
        for k,degree in enumerate(counts):
            choices.append((words[k:]+words[:k],sign*(-1)**(before*(total-before))))
            before+=degree
    least=min(w for w,s in choices)
    signs={s for w,s in choices if w==least}
    return (signs.pop() if len(signs)==1 else 0),least


def variation(row):
    out=Counter()
    for words,c in row.items():
        phase=1
        for i,word in enumerate(words):
            for j,v in enumerate(word):
                for n,u,w in scalar.Q.get(v,()):
                    term=word[:j]+(u,w)+word[j+1:]
                    sg,key=canonical(words[:i]+(term,)+words[i+1:])
                    out[key]+=c*n*phase*sg
                if 3<=v<6:
                    phase=-phase
    return {w:c for w,c in out.items() if c}


@lru_cache(None)
def shapes(length):
    # Three to five adjoint projections, with at least two nontrivial words.
    from itertools import product
    return [s for k in range(3,6) for s in product(range(1,6),repeat=k)
            if sum(s)==length and sum(n>=2 for n in s)>=2]


def draw(length,seed):
    rng=random.Random(seed)
    shape=rng.choice(shapes(length))
    kind="p" if length==11 else rng.choice(("f","pp","pp"))
    slots=list(range(length)); rng.shuffle(slots)
    cuts=np.cumsum((0,)+shape).tolist()
    out=Counter()
    for species,c in scalar.assignments(kind):
        seq=tuple(species[k] for k in slots)
        sg,w=canonical(tuple(seq[a:b] for a,b in zip(cuts,cuts[1:])))
        out[w]+=sg*c
    return {w:c for w,c in out.items() if c}


@lru_cache(1)
def contractions():
    d,_=basis()
    # tr(D_i W D_j) = tr((D_j D_i) W).
    pair=(d[:,None]@d[None,:]).swapaxes(0,1).swapaxes(-1,-2).reshape(52*52,729)
    gram=np.einsum("iab,jba->ij",d,d)
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    return csr_matrix(pair),gram.diagonal().copy()


def values(graphs,f):
    p=f.p
    tensors,gram=contractions()
    assert max(np.abs(tensors).sum(axis=1).A1)*(p-1)<2**63
    inverse=np.array([pow(int(v),-1,p) for v in gram],dtype=np.int64)
    degrees={sum(3<=i<6 for w in g for i in w) for g in graphs}
    assert degrees==set(f.degrees) and degrees<={0,1}
    # Every entry of an unreduced dot product is an integer below 2^53.
    assert 52*(p-1)**2<2**53 and 52**2*(p-1)**2<2**63

    @lru_cache(280)
    def word(w):
        if len(w)==1:
            return (SIGNS[w[0]]*f.L[LETTERS[w[0]]][1])%p
        return word(w[:-1])@word(w[-1:])%p

    @lru_cache(200)
    def projected(w):
        a=word(w)
        vals=(tensors@a.reshape(-1,729).T).T.reshape(*a.shape[:2],52,52)%p
        return (vals*inverse[:,None]%p).astype(np.float64)

    @lru_cache(200)
    def half(words):
        if len(words)==1:
            return projected(words[0])
        return ((half(words[:-1])@projected(words[-1])).astype(np.int64)%p).astype(np.float64)

    result=np.empty((len(graphs),f.ncols),dtype=np.int64)
    for j,g in enumerate(graphs):
        k=len(g)//2
        # Each exact floating-point multiplication is converted separately before summing.
        terms=half(g[:k])*half(g[k:]).swapaxes(-1,-2)
        result[j]=terms.astype(np.int64).sum(axis=(-1,-2)).ravel()%p
    # Compare compressed entries against a direct triple matrix product.
    d,_=basis()
    w=graphs[0][0]
    for i,j in ((0,0),(4,17),(31,14)):
        direct=(d[i]@word(w)%p@d[j]).trace(axis1=-1,axis2=-2)%p
        np.testing.assert_array_equal(projected(w)[:,:,i,j],direct*inverse[i]%p)
    for j in sorted({0,len(graphs)//2,len(graphs)-1}):
        blocks=graphs[j]
        full=projected(blocks[0]).astype(np.int64)
        for block in blocks[1:]:
            full=full@projected(block).astype(np.int64)%p
        np.testing.assert_array_equal(full.trace(axis1=-2,axis2=-1).ravel()%p,result[j])
    half.cache_clear()
    projected.cache_clear()
    word.cache_clear()
    return result


def run(primary,joined,work,workers,draws):
    work.mkdir(parents=True,exist_ok=True)
    old,extra=read(primary),read(joined)
    assert extra["input_sha256"]==hashlib.sha256(primary.read_bytes()).hexdigest()
    p=old["prime"]
    bmn_convention()
    jordan=native(work,p,"jordan")
    closed={27:m.library(work,27,p)}
    result=[]
    for rec,base in zip(old["records"],extra["records"]):
        length=rec["length"]
        mat,labels,alphabet=inherited(rec,base,work,workers,jordan,closed,p)
        frames=[frame(alphabet,v["points"],v["k"],v["degree"],p,v["seed"]) for v in rec["frames"]]
        limit=61 if length==11 else 141-result[0]["minor"]["rank"]
        trace=[]
        first=7501101+4001*length
        for start in range(first,first+draws,12):
            additions=[];seeds=[]
            for seed in range(start,min(start+12,first+draws)):
                expr=variation(draw(length,seed))
                if expr:
                    assert not variation(expr)
                    additions.append(expr);seeds.append(seed)
            if not additions:
                continue
            graphs=sorted({w for row in additions for w in row})
            ids={w:j for j,w in enumerate(graphs)}
            co=coefficient_matrix([[(ids[w],n) for w,n in row.items()] for row in additions],len(graphs),p)
            def evaluate(f):
                return np.asarray(co@values(graphs,f))%p
            with ThreadPoolExecutor(min(len(frames),workers)) as pool:
                data=np.column_stack(list(pool.map(evaluate,frames)))
            mat=np.vstack([mat,data]);labels += [("compression",s) for s in seeds]
            piv=m.minor(mat,p)["rows"];mat=mat[piv];labels=[labels[i] for i in piv]
            assert len(piv)<=limit
            trace.append([seeds[-1],len(graphs),len(piv)])
            print(length,trace[-1],flush=True)
            if len(piv)==limit:
                break
        result.append(dict(length=length,rows=labels,minor=m.minor(mat,p),history=trace))
        m.write(work/"adjoint_compression.json",dict(prime=p,weight=WEIGHT,records=result,
                primary_sha256=hashlib.sha256(primary.read_bytes()).hexdigest(),
                joined_sha256=hashlib.sha256(joined.read_bytes()).hexdigest(),
                classical_cohomology_upper=141-sum(v["minor"]["rank"] for v in result) if len(result)==2 else None))


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--primary",type=Path,default=Path(__file__).with_name("f4")/"primary.json")
    a.add_argument("--joined",type=Path,default=Path(__file__).with_name("f4")/"joined_images.json")
    a.add_argument("--work",type=Path,required=True)
    a.add_argument("--draws",type=int,default=120)
    a.add_argument("--workers",type=int,default=12)
    args=a.parse_args()
    assert __debug__ and min(args.draws,args.workers)>0
    run(**vars(args))
