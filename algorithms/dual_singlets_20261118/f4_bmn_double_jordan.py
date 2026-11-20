"""Classical BMN images from traces with two traceless Jordan projections in F4."""
import argparse
import hashlib
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
import scalar_neighbors as scalar
from f4_bmn_joined_images import recipe, reduce_graphs
from f4_bmn_primary import LETTERS, SIGNS, WEIGHT, bmn_convention
from f4_bmn_trace_images import baseline
from f4_colour import frame, jordan_colours
from f4_dressed_jordan import short_arm
from f4_formal_cycles import coefficient_matrix
from f4_joined_boundaries import native, open_word, vary
from f4_scalar_correction import evaluate_rows
from f4_singlet_characters import singlets
from fortuity.full.letters import Q_letter, is_ferm
from mixed33_exact import read


@lru_cache(150000)
def canonical(parts):
    a,c,b,d = parts
    u,a = short_arm(a)
    v,b = short_arm(b)
    x,a = open_word(a,"jordan")
    y,b = open_word(b,"jordan")
    s = u*v*x*y
    original,rotated = (a,c,b,d),(b,d,a,c)
    rotation = (-1)**(m.nferm((a,c))*m.nferm((b,d)))
    if original==rotated and rotation<0:
        s = 0
    if rotated<original:
        return s*rotation,rotated
    return s,original


def draw(length,seed):
    rng = random.Random(seed)
    kind = "p" if length==11 else rng.choice(("f","pp","pp"))
    slots = list(range(length))
    rng.shuffle(slots)
    sizes = [(a,c,b,length-a-b-c) for a in range(4,length-3) for b in range(4,length-a+1)
             for c in range(length-a-b+1)]
    a,c,b,d = rng.choice(sizes)
    out = Counter()
    for species,k in scalar.assignments(kind):
        for i in species:
            k *= SIGNS[i]
        w = tuple(LETTERS[species[i]] for i in slots)
        s,g = canonical((w[:a],w[a:a+c],w[a+c:a+c+b],w[a+c+b:]))
        if s:
            out[g] += s*k
    return {g:n for g,n in out.items() if n}


def variation(expr, axis):
    out = Counter()
    for parts,n in expr.items():
        phase = 1
        for block,word in enumerate(parts):
            for slot,letter in enumerate(word):
                i,j,f = letter
                if axis==0:
                    terms = [(phase*c,(u,v)) for c,u,v in Q_letter(letter)]
                elif axis in f and axis-1 not in f:
                    g = tuple(sorted(axis-1 if k==axis else k for k in f))
                    terms = [(1,((i,j,g),))]
                else:
                    terms = []
                for c,insert in terms:
                    changed = word[:slot]+insert+word[slot+1:]
                    s,w = canonical(parts[:block]+(changed,)+parts[block+1:])
                    if s:
                        out[w] += n*c*s
                if is_ferm(letter):
                    phase = -phase
    return {g:n for g,n in out.items() if n}


def values(graphs, f):
    # Images here have exterior degree zero or one. At most one factor in a
    # product carries an exterior generator, so ordinary broadcasting suffices.
    degree, = f.degrees
    assert degree in (0,1) and all(m.nferm(g)==degree for g in graphs)
    insertions = jordan_colours()
    flat = insertions.reshape(26,-1)
    transposed = insertions.swapaxes(-1,-2).reshape(26,-1)
    inv = np.array([pow(int(n),-1,f.p) for n in np.einsum("iab,iba->i",insertions,insertions)])
    @lru_cache(None)
    def word(w):
        if not w:
            return np.eye(27,dtype=np.int64)[None,None,:,:]
        return (word(w[:-1])@f.L[w[-1]][1])%f.p
    @lru_cache(None)
    def projected(w):
        a = word(w)
        coordinates = (a.reshape(*a.shape[:2],-1)@transposed.T)%f.p
        return (((coordinates*inv)%f.p)@flat%f.p).reshape(a.shape)
    @lru_cache(None)
    def half(a,c):
        return projected(a) if not c else projected(a)@word(c)%f.p
    assert 27**2*(f.p-1)**2<2**63
    out = np.empty((len(graphs),f.ncols),np.int64)
    for i,(a,c,b,d) in enumerate(graphs):
        end = half(b,d).swapaxes(-1,-2)
        out[i] = np.sum(half(a,c)*end,axis=(-1,-2)).ravel()%f.p
    # Independent direct trace of the matrix product in its specified order.
    for i in sorted({0,len(graphs)//2,len(graphs)-1}):
        a,c,b,d = graphs[i]
        actual = projected(a)@word(c)%f.p@projected(b)%f.p@word(d)%f.p
        np.testing.assert_array_equal(actual.trace(axis1=-2,axis2=-1).ravel()%f.p,out[i])
    return out


def inherited(prior,extra,work,workers,lib,closed,p):
    length = prior["length"]
    dims = singlets("F4",WEIGHT,length+1)
    mat,labels = baseline(prior,dims,p,closed,lib,work,workers,715180+length*73)
    alphabet = set(LETTERS[:3 if length==11 else 6])
    table = {tuple(label):row for label,row in zip(labels,mat)}
    for label in extra["rows"]:
        if label[0]=="previous":
            continue
        channel,spectator,seed = label
        assert channel=="jordan" and spectator
        img = reduce_graphs(vary(recipe(length,channel,spectator,seed),0,channel),channel)
        ev,_ = evaluate_rows([img],alphabet,prior["frames"],p,work,workers,lib,spectator,channel)
        table[tuple(label)] = ev[0]
    matrix = np.stack([table[tuple(label)] for label in extra["rows"]])
    assert int(nmod_mat(matrix[:,extra["minor"]["columns"]].tolist(),p).det())==extra["minor"]["determinant"]
    return matrix, list(extra["rows"]), alphabet


def run(primary,joined,work,workers,draws):
    old, previous = read(primary),read(joined)
    bmn_convention()
    p = old["prime"]
    assert previous["input_sha256"]==hashlib.sha256(primary.read_bytes()).hexdigest()
    work.mkdir(exist_ok=True,parents=True)
    lib = native(work,p,"jordan")
    closed = {27:m.library(work,27,p)}
    records = []
    for rec,extra in zip(old["records"],previous["records"]):
        length = rec["length"]
        assert extra["length"]==length
        mat,labels,alphabet = inherited(rec,extra,work,workers,lib,closed,p)
        frames = [frame(alphabet,v["points"],v["k"],v["degree"],p,v["seed"]) for v in rec["frames"]]
        limit = 61 if length==11 else 141-records[0]["minor"]["rank"]
        history = []
        start = 215871+length*197
        for off in range(0,draws,48):
            batch = range(start+off,start+min(off+48,draws))
            expressions,indices = [],[]
            for seed in batch:
                row=draw(length,seed)
                if not row:
                    continue
                assert all(not variation(row,h) for h in (1,2))
                q=variation(row,0)
                if q:
                    assert not variation(q,0)
                    expressions.append(q)
                    indices.append(seed)
            if not expressions:
                continue
            graphs=sorted({g for r in expressions for g in r})
            lookup={g:i for i,g in enumerate(graphs)}
            coefficients=coefficient_matrix([[(lookup[g],c) for g,c in r.items()] for r in expressions],len(graphs),p)
            def block(f):
                return np.asarray(coefficients@values(graphs,f))%p
            with ThreadPoolExecutor(min(workers,len(frames))) as pool:
                rows = np.column_stack(list(pool.map(block,frames)))
            mat=np.vstack((mat,rows))
            labels += [("two_jordan",i) for i in indices]
            piv=m.minor(mat,p)["rows"]
            mat=mat[piv];labels=[labels[i] for i in piv]
            history.append([batch.stop-1,len(piv)])
            print(length,len(expressions),len(graphs),"rank",len(piv),flush=True)
            assert len(piv)<=limit
            if len(piv)==limit:
                break
        records.append(dict(length=length,rows=labels,minor=m.minor(mat,p),history=history))
        m.write(work/"double_jordan.json",dict(prime=p,primary_sha256=hashlib.sha256(primary.read_bytes()).hexdigest(),
                joined_sha256=hashlib.sha256(joined.read_bytes()).hexdigest(),weight=WEIGHT,records=records,
                classical_cohomology_upper=141-sum(r["minor"]["rank"] for r in records) if len(records)==2 else None))
    return [r["minor"]["rank"] for r in records]


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--primary",type=Path,default=Path(__file__).with_name("f4")/"primary.json")
    a.add_argument("--joined",type=Path,default=Path(__file__).with_name("f4")/"joined_images.json")
    a.add_argument("--work",type=Path,required=True)
    a.add_argument("--workers",type=int,default=8)
    a.add_argument("--draws",type=int,default=480)
    args=a.parse_args()
    assert __debug__ and min(args.workers,args.draws)>0
    print(run(**vars(args)),flush=True)
