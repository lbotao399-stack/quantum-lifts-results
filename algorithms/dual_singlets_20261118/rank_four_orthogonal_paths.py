"""Exhaust the order-ten orthogonal determinants at terminal weight (3,3,2,2,2)."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import hashlib
import json
import math
import subprocess

from flint import fmpz_mat

import mixed_incoming as m
from orthogonal_relations import pack
from rank_four_relations import batch
from rank_four_terminal import bounds, quotient


WEIGHT = (3, 3, 2, 2, 2)


def allowed_words(k):
    return [v for v in combinations_with_replacement(range(5), k)
            if all(v.count(i) <= 1 for i in (0,1))]


def decorations():
    """Two fixed metric edges give one triple, two doubles, or a closed double."""
    for lengths, spectator in (((3,), False), ((2,2), False), ((2,), True)):
        for grouped in product(*(allowed_words(k) for k in lengths)):
            if tuple(sorted(grouped)) != grouped:
                continue
            used = Counter(i for word in grouped for i in word)
            remaining = tuple(n-used[i] for i,n in enumerate(WEIGHT))
            if min(remaining) < 0:
                continue
            paths = list((i,) for i,n in enumerate(remaining) for _ in range(n))
            closed = grouped[0] if spectator else ()
            if not spectator:
                paths.extend(grouped)
            assert len(paths) == 10
            yield tuple(sorted(Counter(paths).items())), closed


def side_choices(word, multiplicity):
    # LL, LR, RR denote the two ends in the determinant's ordered index sets.
    for ll in range(multiplicity+1):
        for rr in range(multiplicity-ll+1):
            lr = multiplicity-ll-rr
            if len(word)%2 == 0 and (ll or rr):
                continue
            if sum(i < 2 for i in word)%2 and max(ll,lr,rr)>1:
                continue
            yield ll,lr,rr


def cases():
    for paths, closed in decorations():
        for sides in product(*(side_choices(word,n) for word,n in paths)):
            if sum(s[0]-s[2] for s in sides):
                continue
            # Interchanging the two index sets transposes the determinant.
            if sides > tuple((rr,lr,ll) for ll,lr,rr in sides):
                continue
            species, left, right, fixed = [], [], [], []
            for (word,_), counts in zip(paths,sides):
                for endpoints, count in zip(((left,left),(left,right),(right,right)),counts):
                    for _ in range(count):
                        start = 2*len(species)
                        species.extend(word)
                        end = 2*len(species)-1
                        endpoints[0].append(start)
                        endpoints[1].append(end)
                        fixed.extend(range(start+1,end))
            if closed:
                start = 2*len(species)
                species.extend(closed)
                fixed.extend((start+1,start+2,start+3,start))
            assert len(left) == len(right) == 10 and sorted(left+right+fixed) == list(range(24))
            yield species, left+right+fixed


def run(work, workers):
    work = work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    case = bounds(3)
    basis, cov = quotient(WEIGHT,case["formal"])
    columns = [[int(v) for v in row] for row in cov.tolist()]
    tab = work/"basis.txt"
    with tab.open("w") as f:
        print(len(basis),file=f)
        for w in basis:
            print(len(w),*(pack(g) for g in w),file=f)
    exe = work/"brauer"
    subprocess.run(["c++","-O3","-std=c++17",str(Path(__file__).with_name("rank_four_brauer.cpp")),
                    "-o",str(exe)],check=True)
    all_cases = list(cases())
    print("determinant templates",len(all_cases),flush=True)
    vectors, witnesses = [], []
    digest = hashlib.sha256()
    nonzero = 0
    for offset in range(0,len(all_cases),workers*8):
        current = all_cases[offset:offset+workers*8]
        with ThreadPoolExecutor(workers) as pool:
            blocks = pool.map(lambda part:batch(exe,tab,"SO",part),
                              [current[j:j+4] for j in range(0,len(current),4)])
            for index,row in enumerate(r for block in blocks for r in block):
                projection = [sum(c*columns[i][j] for i,c in row.items()) for j in range(cov.ncols())]
                digest.update((json.dumps(projection,separators=(",",":"))+"\n").encode())
                if not any(projection):
                    continue
                nonzero += 1
                if fmpz_mat(vectors+[projection]).rank() > len(vectors):
                    vectors.append(projection)
                    sp,slots = current[index]
                    witnesses.append(dict(index=offset+index,species=sp,slots=slots,
                                          raw=sorted(row.items()),singlet=projection))
        print(offset+len(current),"of",len(all_cases),"rank",len(vectors),flush=True)
    assert vectors and len(vectors)<=case["formal"]-case["cartan"]
    minor = next(dict(columns=cols,determinant=det)
                 for cols in combinations(range(cov.ncols()),len(vectors))
                 for det in [int(fmpz_mat([[v[i] for i in cols] for v in vectors]).det())] if det)
    m.write(work/"orthogonal_paths.json",
            dict(weight=WEIGHT,group="SO(9)",level=30,field="Q",formal=case["formal"],cartan=case["cartan"],
                 templates=len(all_cases),terms_per_template=math.factorial(10),nonzero_projections=nonzero,
                 projection_sha256=digest.hexdigest(),relations=witnesses,rank=len(vectors),
                 minor=minor,classical=case["formal"]-len(vectors)))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    assert __debug__
    run(**vars(p.parse_args()))
