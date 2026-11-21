"""Lift four five-scalar trace families to polynomial first corrections on arbitrary SU(3) matrices."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from su3_level32_layers import content, formal, joint
from su3_level32_levi import extensions
from su3_level32_scalar_descent import physical, restored, scalar_splits, serialized
from verify import add_scaled
import content_parents


FAMILIES=("040404060e16","040606060c14","040406060c16","040406060e14")
PAIRS=extensions()
BARE={marked:k for k,(marked,_) in PAIRS.items()}


def scalar_split(word):
    count=content(word).count(0)
    return {w:c for w,c in formal.q(word).items() if content(w).count(0)==count+1}


def normalize(row):
    den=lcm(*(c.denominator for c in row.values()))
    nums={w:int(den*c) for w,c in row.items()}
    scale=gcd(*nums.values())
    assert scale
    return {w:c//scale for w,c in nums.items()},Fraction(scale,den)


def factor_rows(original):
    by_content=defaultdict(dict)
    for w,c in original.items():
        by_content[content(w)][w]=c
    cores=defaultdict(dict)
    for row in by_content.values():
        common=Counter(next(iter(row)))
        for w in row:
            common &= Counter(w)
        traces=tuple(sorted(common.elements()))
        assert len(traces)==2 and bytes(2) in traces and not any(scalar_split((t,)) for t in traces)
        for w,c in row.items():
            rest=list(w)
            for t in traces:
                rest.remove(t)
            sign,canon=formal.multi(traces+tuple(rest))
            assert canon==w and sign
            add_scaled(cores[traces],{tuple(rest):sign},c)
    return cores


def cofactor_parents(core):
    bare_words={bytes(sorted(BARE.get(k,k) for k in content(w) if k)) for w in core}
    assert len(bare_words)==1
    bare=bare_words.pop()
    assert len(bare)==5 and all(content(w).count(0)==2 and sum(k in BARE for k in content(w))==1 for w in core)
    marked={bytes(sorted([0]+[PAIRS[k][0] if i in choices else k for i,k in enumerate(bare)]))
            for choices in combinations(range(len(bare)),2)}
    return sorted({w for t in marked for w in content_parents.all_words(t)})


def solve_core(core,engine,work):
    rhs,scale=normalize(core)
    parents=cofactor_parents(core)
    images=[scalar_split(w) for w in parents]
    words=sorted(set(rhs).union(*(image.keys() for image in images)))
    keys={}
    polys={}
    for pos in range(0,len(words),40):
        rows=words[pos:pos+40]
        # A common tr(S^4) puts every word in the eleven-letter exporter.
        # This even nonzero scalar polynomial acts injectively on the exterior coefficients.
        padded=[tuple(sorted((bytes(4),)+w)) for w in rows]
        for w,v in zip(rows,joint.E.poly_iter(padded,engine,work),strict=True):
            polys[w]={keys.setdefault(k,len(keys)):n for k,n in v}

    def polynomial(row):
        value=Counter()
        for w,c in row.items():
            for k,n in polys[w].items():
                value[k]+=c*n
        return {k:n for k,n in value.items() if n}

    offsets=[0]
    indices=[]
    entries=[]
    for row in images:
        for i,c in sorted(polynomial(row).items()):
            indices.append(i)
            entries.append(c)
        offsets.append(len(indices))
    matrix=csr_matrix((np.array(entries,dtype=np.int64),indices,offsets),shape=(len(parents),len(keys)))
    terms=polynomial(rhs)
    values=[-terms.get(i,0) for i in range(len(keys))]
    info,solution=physical.certify(matrix,values,lambda _:None)
    assert info["kind"]=="exact_solution"
    check=rhs.copy()
    correction={}
    for i,c in solution.items():
        add_scaled(check,images[i],c)
        correction[parents[i]]=scale*c
    assert not polynomial(check)
    stats=dict(rank=info,trace_parents=len(parents),trace_words=len(words),colour_monomials=len(keys),
               nonzero_solution_terms=len(correction),integer_rhs_terms=len(rhs))
    return correction,stats


def solve(folder,engine,work):
    work.mkdir(parents=True,exist_ok=True)
    records={}
    combined={}
    selected_contents=set()
    for bare in FAMILIES:
        source=joint.s.read(folder/(bare+".json.gz"))
        assert selected_contents.isdisjoint(source)
        selected_contents.update(source)
        original={w:c for rows in source.values() for w,c in restored(rows).items()}
        correction={}
        records[bare]=[]
        for factors,core in sorted(factor_rows(original).items()):
            small,stats=solve_core(core,engine,work)
            parity=(-1)**sum(formal.par(t) for t in factors)
            for w,c in small.items():
                sign,canon=formal.multi(factors+w)
                if sign:
                    add_scaled(correction,{canon:sign*parity},c)
            records[bare].append(dict(spectators=[t.hex() for t in factors],**stats))
        check=original.copy()
        for w,c in correction.items():
            assert content(w).count(0)==4
            add_scaled(check,scalar_splits(w),c)
        by_content=defaultdict(dict)
        for w,c in check.items():
            by_content[content(w).hex()][w]=c
        assert set(by_content)<=set(source)
        checks={}
        for t,row in sorted(by_content.items()):
            verified=physical.polynomials([row],engine,work)[0]
            assert not verified["nonzero_monomials"]
            checks[t]=verified
        records[bare]=dict(cores=records[bare],target_contents=sorted(source),parent_terms=len(correction),
                           full_three_colour_checks=checks)
        add_scaled(combined,correction,1)
        print(bare,len(correction),len(source),sum(v["collected_monomials"] for v in checks.values()),flush=True)
    joint.s.write(work/"factored_delta.json.gz",serialized(combined))
    joint.s.write(work/"factored_descent.json",dict(families=records,delta_terms=len(combined),
                   target_contents=len(selected_contents)))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder",type=Path,required=True)
    p.add_argument("--engine",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    solve(args.folder,args.engine.resolve(),args.work)
