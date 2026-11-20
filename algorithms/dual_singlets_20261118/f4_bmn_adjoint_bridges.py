"""BMN boundaries from two traces joined by two independent adjoint contractions in F4."""
import argparse
import hashlib
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np

import mixed_incoming as m
import scalar_neighbors as scalar
from f4_bmn_primary import LETTERS, SIGNS, WEIGHT, bmn_convention
from f4_bmn_double_jordan import inherited
from f4_colour import basis, frame
from f4_formal_cycles import coefficient_matrix
from f4_joined_boundaries import native
from fortuity.full.letters import Q_letter, is_ferm
from mixed33_exact import read


@lru_cache(160000)
def canonical(parts):
    a,b,c,d = parts
    pa,pb,pc,pd = (m.nferm((w,)) for w in parts)
    swap_labels = (-1)**(pa*pb+pc*pd)
    swap_traces = (-1)**((pa+pb)*(pc+pd))
    images = [(parts,1),((b,a,d,c),swap_labels),((c,d,a,b),swap_traces),
              ((d,c,b,a),swap_labels*swap_traces)]
    canonical_word = min(w for w,s in images)
    signs = {s for w,s in images if w==canonical_word}
    return (signs.pop() if len(signs)==1 else 0),canonical_word


def draw(length,seed):
    rng = random.Random(seed)
    kind = "p" if length==11 else rng.choice(("f","pp","pp"))
    slots = list(range(length))
    rng.shuffle(slots)
    shapes = [(a,b,c,length-a-b-c) for a in (1,2,3) for b in range(a,5) for c in (1,2,3)
              if 4<=a+b<=6 and 4<=length-a-b<=7 and c<=length-a-b-c<=5]
    n,k,j,_ = rng.choice(shapes)
    result = Counter()
    for species,co in scalar.assignments(kind):
        for i in species:
            co *= SIGNS[i]
        w = tuple(LETTERS[species[i]] for i in slots)
        s,g = canonical((w[:n],w[n:n+k],w[n+k:n+k+j],w[n+k+j:]))
        if s:
            result[g] += s*co
    return {g:v for g,v in result.items() if v}


def variation(expr,axis):
    result = Counter()
    for parts,coef in expr.items():
        phase=1
        for j,w in enumerate(parts):
            for i,letter in enumerate(w):
                a,b,f=letter
                if axis==0:
                    terms=[(phase*c,(x,y)) for c,x,y in Q_letter(letter)]
                elif axis in f and axis-1 not in f:
                    terms=[(1,((a,b,tuple(sorted(axis-1 if k==axis else k for k in f))),))]
                else:
                    terms=[]
                for c,letters in terms:
                    s,g=canonical(parts[:j]+(w[:i]+letters+w[i+1:],)+parts[j+1:])
                    if s:
                        result[g]+=s*c*coef
                if is_ferm(letter):
                    phase=-phase
    return {g:v for g,v in result.items() if v}


def values(graphs,f):
    degree,=f.degrees
    assert degree in (0,1) and all(m.nferm(g)==degree for g in graphs)
    generators,_=basis()
    gram=np.einsum("iab,jba->ij",generators,generators)
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    dual=np.array([pow(int(n),-1,f.p) for n in gram.diagonal()],np.int64)
    weights=dual[:,None]*dual[None,:]%f.p
    assert 729*(f.p-1)**2<2**53 and 52**2*(f.p-1)**2<2**63

    @lru_cache(None)
    def word(w):
        if len(w)==1:
            return f.L[w[0]][1]
        return word(w[:-1])@f.L[w[-1]][1]%f.p

    @lru_cache(None)
    def factors(w,transpose):
        mat=generators[None,None,:,:,:]@word(w)[:,:,None,:,:]%f.p
        if transpose:
            mat=mat.swapaxes(-1,-2)
        return np.ascontiguousarray(mat.reshape(*mat.shape[:3],729),dtype=np.float64)

    @lru_cache(None)
    def joined(a,b):
        # Each dot product adds 729 nonnegative integers. The asserted 2^53
        # bound makes every binary64 multiply and addition exact.
        val=(factors(a,False)@factors(b,True).swapaxes(-1,-2)).astype(np.int64)%f.p
        return val.astype(np.uint32),(val*weights%f.p).astype(np.uint32)

    out=np.empty((len(graphs),f.ncols),np.int64)
    for i,(a,b,c,d) in enumerate(graphs):
        left=joined(a,b)[1]
        right=joined(c,d)[0]
        out[i]=np.sum(np.multiply(left,right,dtype=np.int64),axis=(-1,-2)).ravel()%f.p
    a,b,c,d=graphs[0]
    u,v=0,17
    for first,second in [(a,b),(c,d)]:
        control=generators[u]@word(first)%f.p@generators[v]%f.p@word(second)%f.p
        np.testing.assert_array_equal(control.trace(axis1=-2,axis2=-1)%f.p,joined(first,second)[0][:,:,u,v])
    return out


def run(primary,joined,work,workers,draws):
    old,prev=read(primary),read(joined)
    p=old["prime"]
    bmn_convention()
    assert prev["input_sha256"]==hashlib.sha256(primary.read_bytes()).hexdigest()
    work.mkdir(exist_ok=True,parents=True)
    jordan=native(work,p,"jordan")
    closed={27:m.library(work,27,p)}
    results=[]
    for rec,base in zip(old["records"],prev["records"]):
        length=rec["length"]
        assert length==base["length"]
        mat,labels,alph=inherited(rec,base,work,workers,jordan,closed,p)
        frames=[frame(alph,v["points"],v["k"],v["degree"],p,v["seed"]) for v in rec["frames"]]
        limit=61 if length==11 else 141-results[0]["minor"]["rank"]
        first=6178501+1531*length
        history=[]
        for offset in range(0,draws,24):
            seeds=range(first+offset,first+min(offset+24,draws))
            images,extras=[],[]
            for seed in seeds:
                row=draw(length,seed)
                assert all(not variation(row,k) for k in (1,2))
                img=variation(row,0)
                if img:
                    assert not variation(img,0)
                    images.append(img)
                    extras.append(("adjoint_bridges",seed))
            if not images:
                continue
            graphs=sorted({g for r in images for g in r})
            ids={g:i for i,g in enumerate(graphs)}
            co=coefficient_matrix([[(ids[g],c) for g,c in r.items()] for r in images],len(graphs),p)
            def evaluate(f):
                return np.asarray(co@values(graphs,f))%p
            with ThreadPoolExecutor(min(6,workers,len(frames))) as pool:
                additions=np.column_stack(list(pool.map(evaluate,frames)))
            mat=np.vstack((mat,additions))
            labels+=extras
            piv=m.minor(mat,p)["rows"]
            mat=mat[piv];labels=[labels[i] for i in piv]
            assert len(piv)<=limit
            history.append([seeds.stop-1,len(piv)])
            print(length,len(graphs),"rank",len(piv),flush=True)
            if len(piv)==limit:
                break
        results.append(dict(length=length,rows=labels,minor=m.minor(mat,p),history=history))
        m.write(work/"adjoint_bridges.json",dict(prime=p,weight=WEIGHT,
            primary_sha256=hashlib.sha256(primary.read_bytes()).hexdigest(),
            joined_sha256=hashlib.sha256(joined.read_bytes()).hexdigest(),records=results,
            classical_cohomology_upper=141-sum(r["minor"]["rank"] for r in results) if len(results)==2 else None))
    return [r["minor"]["rank"] for r in results]


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--primary",type=Path,default=Path(__file__).with_name("f4")/"primary.json")
    a.add_argument("--joined",type=Path,default=Path(__file__).with_name("f4")/"joined_images.json")
    a.add_argument("--work",type=Path,required=True)
    a.add_argument("--workers",type=int,default=12)
    a.add_argument("--draws",type=int,default=144)
    args=a.parse_args()
    assert __debug__ and min(args.workers,args.draws)>0
    print(run(**vars(args)),flush=True)
