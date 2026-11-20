"""Classical bounds for the bilinear multiples of possible F4 second- and third-order sources."""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np

import mixed_incoming as m
from f4_formal_cycles import coefficient_matrix, formal_columns, probe
from f4_singlet_characters import singlets
from mixed33_exact import exact, read
from mixed_formal import degree
from mixed_middle import build, pin


def rational_rank(images, dictionary):
    blocks = defaultdict(list)
    for row in images:
        if row:
            labels = {degree(dictionary[i]) for i,c in row}
            assert len(labels)==1
            blocks[labels.pop()].append(dict(row))
    return sum(len(exact(part)[0]) for part in blocks.values())


def minor(x,p):
    return m.minor(x,p) if x.size and np.any(x) else dict(rank=0)


def run(work, workers, cache=None, small=False):
    work.mkdir(parents=True,exist_ok=True)
    p = 1000003
    lib = m.library(work,27,p)
    records = []
    for order,t in ((2,3),(3,2)):
        for factors in reversed(range(1,t+1)):
            if small and factors==1: continue
            spin = t-factors
            power = 5-factors
            weight = (spin,spin)+(5-t,)*3
            n = sum(weight)-1
            dims = singlets("F4",weight,n)
            ambient = sum(dims.values())
            ends = []
            for length in (n-1,n):
                parent = work/f"p{power}-t{spin}-n{length}"
                parent.mkdir(exist_ok=True)
                path = parent/"plan.json.gz"
                if cache and power==4 and length==n-1:
                    candidate = cache/"dual_mixed_middle_v1_20261120"/f"t{spin}"/"plan.json.gz"
                    if candidate.is_file(): path=candidate
                if not path.is_file():
                    path=build(spin,parent,workers,length=length,power=power)
                z=read(path)
                assert tuple(z["mu"])==weight and z["length"]==length
                images=z["images"]
                words=list(map(m.decode,z["target_dictionary"]))
                d=singlets("F4",weight,length+1)
                fd={m.nferm(w) for w in words}
                if words:
                    values,_=probe(images,words,{k:d[k] for k in fd},p,lib,workers,1507032+length*229)
                    rank=minor(values,p)
                else: rank=dict(rank=0)
                ends.append(dict(plan_sha256=pin(path),rank=rank))
            sw=list(map(m.decode,z["source_dictionary"]))
            src,_=probe(z["sources"],sw,dims,p,lib,workers,1250622+n*229)
            formal_rank=rational_rank(images,words)
            cols=formal_columns(images,p)
            assert len(cols)==formal_rank
            formal=coefficient_matrix(images,len(words),p)[:,cols].toarray()
            joint=minor(np.column_stack((src,formal)),p)
            closed=joint["rank"]-formal_rank
            total=ambient-ends[0]["rank"]["rank"]-ends[1]["rank"]["rank"]
            remaining=ambient-ends[1]["rank"]["rank"]-closed
            assert total>=0 and remaining>=0
            records.append(dict(order=order,S_factors=factors,weight=weight,base_length=n,cochain=ambient,
                                incoming=ends[0],outgoing=ends[1],classical_upper=total,formal_rank=formal_rank,
                                restricted_formal_closed_lower=closed,joint=joint,
                                kernel_modulo_formal_upper=remaining))
            m.write(work/"dressed_sources.json",dict(prime=p,records=records))
            print("order",order,"S factors",factors,"cohomology upper",total,
                  "kernel modulo formal upper",remaining,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--cache",type=Path)
    p.add_argument("--small",action="store_true")
    assert __debug__
    run(**vars(p.parse_args()))
