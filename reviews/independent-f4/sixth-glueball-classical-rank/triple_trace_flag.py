"""Append three-trace classical images with zero-score primitives to the F4 flag minor."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np

from double_trace_flag import Evaluation, necklace
from flag_boundaries import HERE, flags
from single_trace_flag import extend, native26, trace_values


SHAPES = [(a,b,14-a-b) for a in range(5,11) for b in range(2,a+1)
          if 2 <= 14-a-b <= b]


def terms(shape,seed):
    digits,signs = flags()
    rng = np.random.Generator(np.random.PCG64(3604433+1009*seed+100*shape[0]+shape[1]))
    groups = np.split(rng.permutation(14),np.cumsum(shape)[:-1])
    packed = np.zeros(len(signs),dtype=np.int64)
    phase = signs.copy()
    for n,slots in zip(shape,groups):
        word = digits[slots[0]].copy()
        for i in slots[1:]:
            word = 4*word+digits[i]
        key,s = necklace(word,n)
        packed = (packed << (2*n)) + key
        phase *= s
    keys,inverse = np.unique(packed,return_inverse=True)
    coeff = np.bincount(inverse,weights=phase).astype(np.int64)
    return keys[coeff!=0],coeff[coeff!=0]


def direct_value(shape,seed,number,prime,rep):
    """Integer evaluation of the uncollected flag, keeping oriented trace words."""
    digits,signs = flags()
    permutation = np.random.Generator(np.random.PCG64(3604433+1009*seed+100*shape[0]+shape[1])).permutation(14)
    coeff = np.random.Generator(np.random.PCG64(144336+number)).integers(-3,4,size=(4,52))
    x = (coeff @ rep.reshape(52,-1)).reshape(4,26,26)%prime
    cache = {(): np.eye(26,dtype=np.int64)}
    def matrix(w):
        if w not in cache:
            cache[w] = matrix(w[:-1])@x[w[-1]]%prime
        return cache[w]
    value = signs.astype(np.int64)
    for slots in np.split(permutation,np.cumsum(shape)[:-1]):
        words,where = np.unique(digits[slots].T,axis=0,return_inverse=True)
        half = len(slots)//2
        traces = np.array([np.sum(matrix(tuple(w[:half]))*matrix(tuple(w[half:])).T)%prime for w in words])
        value = value*traces[where]%prime
    return int(value.sum())%prime


class Triple:
    def __init__(self,rows):
        self.blocks = []
        all_words = {n:[] for n in range(2,11)}
        for shape in SHAPES:
            name = "T"+".".join(map(str,shape))
            chosen = [(i,seed) for i,(k,seed) in enumerate(rows) if k==name]
            if not chosen:
                continue
            polys = [terms(shape,seed) for _,seed in chosen]
            keys = np.unique(np.concatenate([a for a,_ in polys]))
            ids = np.concatenate([np.searchsorted(keys,a) for a,_ in polys])
            coeff = np.concatenate([v for _,v in polys])
            ends = np.cumsum([0]+[len(a) for a,_ in polys])
            active = np.flatnonzero(np.diff(ends))
            arrays = []
            for n in reversed(shape):
                part,where = np.unique(keys & ((1<<(2*n))-1),return_inverse=True)
                arrays.append((n,part,where))
                all_words[n].append(part)
                keys = keys >> (2*n)
            self.blocks.append((np.array([i for i,_ in chosen])[active],ends[active],arrays,ids,coeff))
        self.words = {n:np.unique(np.concatenate(words)) for n,words in all_words.items() if words}
        for _,_,arrays,_,_ in self.blocks:
            for i,(n,part,where) in enumerate(arrays):
                arrays[i] = (n,np.searchsorted(self.words[n],part)[where])
        self.count = len(rows)
        self.representation = native26()

    def point(self,number,prime,integer=False):
        coefficients = np.random.Generator(np.random.PCG64(144336+number)).integers(-3,4,size=(4,52))
        x = (coefficients @ self.representation.reshape(52,-1)).reshape(4,26,26)%prime
        if not integer:
            x = x.astype(np.float64)
        tables = {1:x}
        for n in range(2,(max(self.words)+1)//2+1):
            with np.errstate(divide="ignore",over="ignore",invalid="ignore"):
                raw = tables[n-1][:,None] @ x[None]
            assert 0 <= raw.min() and raw.max() <= 26*(prime-1)**2 < 2**53
            tables[n] = (raw%prime).reshape(-1,26,26)
        values = {n:trace_values(tables,words,n,prime) for n,words in self.words.items()}
        out = np.zeros(self.count,dtype=np.int64)
        for selected,starts,parts,indices,coeff in self.blocks:
            if len(selected):
                v = np.ones(len(parts[0][1]),dtype=np.int64)
                for n,pos in parts:
                    v = v*values[n][pos]%prime
                out[selected] = np.add.reduceat(coeff*v[indices],starts)%prime
        return out


def calculate(args):
    reference = json.loads(args.replay.read_bytes()) if args.replay else None
    if reference:
        recipes = reference["rows"]
        points = reference["points"]
    else:
        recipes = json.loads((HERE/"double_selected.json").read_bytes())["rows"]
        recipes += [["T"+".".join(map(str,shape)),s] for shape in SHAPES for s in range(args.seeds)]
        points = list(range(704))
    old = [r for r in recipes if not r[0].startswith("T")]
    new = [r for r in recipes if r[0].startswith("T")]
    prior = Evaluation(old)
    triples = Triple(new)
    labels = old+new
    assert sorted(map(tuple,recipes)) == sorted(map(tuple,labels))
    old_first = {tuple(r):i for i,r in enumerate(labels)}
    reorder = [old_first[tuple(r)] for r in recipes]
    control = triples.point(points[0],args.prime,True)
    assert np.array_equal(triples.point(points[0],args.prime),control)
    j = int(np.flatnonzero(control)[0])
    group,seed = new[j]
    shape = tuple(map(int,group[1:].split(".")))
    assert direct_value(shape,seed,points[0],args.prime,triples.representation) == control[j]

    def point(p):
        return np.concatenate((prior.point(p,args.prime),triples.point(p,args.prime)))[reorder]

    pivots, kept, rows = {}, [], []
    determinant = 1
    with ThreadPoolExecutor(args.workers) as pool:
        for start in range(0,len(points),args.workers):
            batch = points[start:start+args.workers]
            for p,values in zip(batch,pool.map(point,batch)):
                pivot = extend(pivots,values,args.prime)
                if pivot:
                    i,v = pivot
                    determinant = determinant*v % args.prime
                    rows.append(i)
                    kept.append(p)
            if start//args.workers%4==0:
                print("points",start+len(batch),"rank",len(kept),flush=True)
            assert len(kept)<=703
            if len(kept)==703:
                break
    if reference:
        if len(rows)==len(recipes):
            assert kept == points
            inversions = sum(i>j for k,i in enumerate(rows) for j in rows[k+1:])
            determinant = determinant*(-1)**(inversions%2)%args.prime
        else:
            determinant = 0
        return dict(prime=args.prime,rank=len(rows),determinant=determinant,points=points,rows=recipes)
    return dict(prime=args.prime,rank=len(rows),determinant=determinant,points=kept,rows=[recipes[i] for i in rows])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime",type=int,choices=(1009,1013),required=True)
    parser.add_argument("--seeds",type=int,default=128)
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--replay",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    if not __debug__ or min(args.seeds,args.workers)<1:
        parser.error("enable assertions and use positive sizes")
    record = calculate(args)
    args.out.write_text(json.dumps(record,indent=2)+"\n")
    print("rank",record["rank"],"determinant",record["determinant"],flush=True)
