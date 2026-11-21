"""Extend an invariant scalar-family correction in the six-scalar equations modulo I squared."""
import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from math import lcm
from pathlib import Path

from mixed33_exact import exact
from su3_square_zero_normal import SCALARS, disjoint_pairs, filtered_residual, formal


def scalar_count(word):
    return sum(k in SCALARS for t in word for k in t)


@lru_cache(maxsize=20000)
def high_variation(word,minimum_first):
    return {w:n for w,n in formal.q(word).items()
            if scalar_count(w)>=6 and sum(t.count(0) for t in w)>=minimum_first
            and sum(map(disjoint_pairs,w))<2}


def connected_parents(seeds,minimum_first):
    """Explore all degree-six-or-higher images within a fixed multiset of trace weights."""
    frontier = set(seeds)
    seen = set(frontier)
    examined = set()
    columns = {}
    while frontier:
        next_words = set()
        for w in frontier:
            for parent in formal.parents(w):
                if parent in examined:
                    continue
                examined.add(parent)
                if sum(map(disjoint_pairs,parent))!=1:
                    continue
                image = high_variation(parent,minimum_first)
                if w not in image:
                    examined.discard(parent)
                    continue
                columns[parent] = image
                next_words.update(image)
        frontier = next_words-seen
        seen.update(frontier)
    return columns


def correction_block(residual,minimum_first):
    columns = connected_parents(residual,minimum_first)
    support = {w:set(v) for w,v in columns.items()}
    incident = defaultdict(list)
    for parent in sorted(columns):
        for w in columns[parent]:
            incident[w].append(parent)
    queue = deque(p for p in sorted(columns) if len(support[p])==1)
    peeled = []
    while queue:
        parent = queue.popleft()
        if len(support[parent])!=1:
            continue
        w, = support[parent]
        peeled.append((w,parent))
        for other in incident[w]:
            support[other].remove(w)
            if len(support[other])==1:
                queue.append(other)
    removed = {w for w,_ in peeled}
    parents = sorted(p for p,s in support.items() if s)
    words = sorted(({w for p in parents for w in support[p]} | set(residual))-removed)
    solution = {}
    if words and parents and any(w in residual for w in words):
        ids = {w:i for i,w in enumerate(words)}
        den = lcm(*(residual[w].denominator for w in words if w in residual))
        rhs = {w:int(residual[w]*den) for w in words if w in residual}
        arrays = [{ids[w]:columns[parent][w] for w in support[parent]} for parent in parents]
        arrays.append({ids[w]:n for w,n in rhs.items()})
        independent,_,relations = exact(arrays)
        if len(parents) not in independent:
            relation = relations[-1]
            solution = {parents[i]:Fraction(n,den*relation[len(parents)]) for i,n in relation.items()
                        if i<len(parents)}
    check = Counter(residual)
    for parent,c in solution.items():
        for w,n in columns[parent].items():
            check[w] += c*n
    for w,parent in reversed(peeled):
        if check[w]:
            c = -Fraction(check[w],columns[parent][w])
            solution[parent] = c
            for v,n in columns[parent].items():
                check[v] += c*n
    remaining = {w:c for w,c in check.items() if c}
    return solution,remaining,len(columns),len(words)


def run(family, destination, build,minimum_first):
    rows = filtered_residual(family)
    high = {t:v for t,v in rows.items() if sum(k in SCALARS for k in t)>=6 and t.count(0)>=minimum_first}
    if build:
        solution = Counter()
        blocks = defaultdict(dict)
        for r in high.values():
            for w,c in r.items():
                blocks[formal.key(w)][w] = c
        for i,(key,block) in enumerate(sorted(blocks.items(),key=lambda item:(len(item[1]),item[0]))):
            coeff,remainder,parents,core = correction_block(block,minimum_first)
            solution.update(coeff)
            print(i,key,len(block),len(remainder),parents,core,flush=True)
        polynomial = [[formal.display(w),str(c)] for w,c in sorted(solution.items()) if c]
        destination.write_bytes(gzip.compress((json.dumps(polynomial,separators=(",",":"))+"\n").encode(),mtime=0))
    polynomial = json.loads(gzip.decompress(destination.read_bytes()))
    corrected = {t:Counter(row) for t,row in high.items()}
    for raw,c in polynomial:
        w = formal.unpack(raw)
        assert sum(map(disjoint_pairs,w))>=1
        for v,n in formal.q(w).items():
            if scalar_count(v)>=6 and sum(t.count(0) for t in v)>=minimum_first and sum(map(disjoint_pairs,v))<2:
                t = bytes(sorted(b"".join(v)))
                corrected.setdefault(t,Counter())[v] += Fraction(c)*n
    results = [{"content":t.hex(),"before":len(high.get(t,{})),"after":sum(bool(c) for c in v.values())}
               for t,v in sorted(corrected.items())]
    return dict(level=32,minimum_first_flavour=minimum_first,family_sha256=hashlib.sha256(family.read_bytes()).hexdigest(),
                normal_correction_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
                correction_terms=len(polynomial),contents=results,
                solved_initial_contents=sum(r["before"]>0 and r["after"]==0 for r in results),
                initial_terms=sum(r["before"] for r in results),remaining_terms=sum(r["after"] for r in results))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--family",type=Path,required=True)
    p.add_argument("--correction",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--build",action="store_true")
    p.add_argument("--minimum-first-flavour",type=int,default=0)
    a=p.parse_args()
    if not __debug__:p.error("enable assertions")
    result=run(a.family,a.correction,a.build,a.minimum_first_flavour)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
