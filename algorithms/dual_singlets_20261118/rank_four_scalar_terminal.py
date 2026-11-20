"""Orthogonal determinant relations on the fourteen-letter level-30 terminal singlets."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from math import gcd
from pathlib import Path
import random
import subprocess

from flint import fmpz_mat

import mixed_incoming as m
from orthogonal_relations import lower, pack
from rank_four_relations import batch
from rank_four_terminal import bounds
from terminal_commutators import words


def trace_shape(word):
    return tuple(sorted((sum(v),v[0]+v[1]) for v in word))


def split_quotient(weight):
    basis = words(weight)
    idx = {w:i for i,w in enumerate(basis)}
    blocks = defaultdict(list)
    for i,w in enumerate(basis):
        blocks[trace_shape(w)].append(i)
    equations = defaultdict(list)
    for i,j in ((0,1),(2,3),(3,4)):
        mu = list(weight); mu[i] += 1; mu[j] -= 1
        for w in words(tuple(mu)):
            row = {idx[v]:c for v,c in lower(w,i,j).items() if c}
            if row:
                key = trace_shape(w)
                assert all(trace_shape(basis[k])==key for k in row)
                equations[key].append(row)
    coordinates = [[] for _ in basis]
    blocks_used = []
    for key,indices in sorted(blocks.items()):
        rows = equations[key]
        a = fmpz_mat(len(rows),len(indices))
        loc = {v:i for i,v in enumerate(indices)}
        for i,row in enumerate(rows):
            for j,c in row.items():
                a[i,loc[j]] = c
        nullity = a.ncols()-a.rank()
        if not nullity:
            continue
        null,k = a.nullspace()
        assert k==nullity
        for j in range(k):
            col = [int(null[i,j]) for i in range(len(indices))]
            g = gcd(*col)
            vec = fmpz_mat([[c//g] for c in col])
            assert (a*vec).is_zero()
            for v,row in enumerate(coordinates):
                row.append(int(vec[loc[v],0]) if v in loc else 0)
        blocks_used.append(dict(shape=key,weight_words=len(indices),singlets=k))
    return basis,coordinates,blocks_used


def sample(rng, weight):
    species = [i for i,n in enumerate(weight) for _ in range(n)]
    rng.shuffle(species)
    # Ten free-ended paths and closed spectator traces account for four fixed metric edges.
    loops = rng.choice(((),(),(),(2,),(2,),(4,),(2,2)))
    lengths = [1]*10
    for _ in range(4-sum(loops)):
        lengths[rng.randrange(10)] += 1
    paths, fixed, cursor = [], [], 0
    for count in (*lengths,*loops):
        start = 2*cursor
        cursor += count
        end = 2*cursor-1
        fixed.extend(range(start+1,end))
        if len(paths)<10:
            paths.append((count,start,end))
        else:
            fixed.extend((end,start))
    # Even paths on a single side antisymmetrize to zero.
    odd = [i for i,(k,_,_) in enumerate(paths) if k%2]
    rng.shuffle(odd)
    pairs = rng.randrange(len(odd)//2+1)
    ll,rr = set(odd[:pairs]),set(odd[pairs:2*pairs])
    left,right = [],[]
    for i,(_,start,end) in enumerate(paths):
        if i in ll: left.extend((start,end))
        elif i in rr: right.extend((start,end))
        else: left.append(start); right.append(end)
    assert len(left)==len(right)==10 and sorted(left+right+fixed)==list(range(28))
    return species,left+right+fixed


def run(work, workers, attempts):
    work = work.resolve(); work.mkdir(parents=True,exist_ok=True)
    case = bounds(1)
    weight = case["weight"]
    basis,cols,shape = split_quotient(weight)
    assert len(basis)==6702 and len(cols[0])==case["formal"]==18
    table = work/"basis.txt"
    with table.open("w") as f:
        print(len(basis),file=f)
        for w in basis: print(len(w),*(pack(g) for g in w),file=f)
    exe = work/"brauer"
    subprocess.run(["c++","-O3","-std=c++17","-DRANK_FOUR_LETTERS=14",
                    str(Path(__file__).with_name("rank_four_brauer.cpp")),"-o",str(exe)],check=True)
    rng = random.Random(871143)
    rows,selected = [],[]
    examined = 0
    for start in range(0,attempts,workers*16):
        cases = [sample(rng,weight) for _ in range(min(attempts-start,workers*16))]
        with ThreadPoolExecutor(workers) as pool:
            chunks = pool.map(lambda v:batch(exe,table,"SO",v),[cases[i:i+4] for i in range(0,len(cases),4)])
            for i,raw in enumerate(r for c in chunks for r in c):
                examined += 1
                v = [sum(c*cols[j][h] for j,c in raw.items()) for h in range(18)]
                if any(v) and fmpz_mat(rows+[v]).rank()>len(rows):
                    rows.append(v)
                    sp,slots = cases[i]
                    selected.append(dict(trial=start+i,species=sp,slots=slots,raw=sorted(raw.items()),singlet=v))
                    print("rank",len(rows),"trial",start+i,flush=True)
                    if len(rows)==12: break
        minor = {}
        if rows:
            rref = fmpz_mat(rows).rref()[0]
            pivots = [next(j for j in range(rref.ncols()) if rref[i,j]) for i in range(len(rows))]
            det = int(fmpz_mat([[v[j] for j in pivots] for v in rows]).det())
            assert det
            minor = dict(columns=pivots,determinant=det)
        m.write(work/"scalar_terminal.json",dict(weight=weight,formal=18,cartan=6,lower=6,upper=18-len(rows),
                seed=871143,examined=examined,terms_per_relation=3628800,highest_weight_blocks=shape,relations=selected,
                minor=minor))
        if len(rows)==12: break
    print("SO(9) fourteen-letter bounds",6,18-len(rows),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--attempts",type=int,default=12000)
    assert __debug__
    run(**vars(p.parse_args()))
