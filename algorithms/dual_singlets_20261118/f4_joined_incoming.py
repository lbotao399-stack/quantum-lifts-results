"""Add classical boundaries from nine-letter F4 covariants to the universal closed span."""
import argparse
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from f4_colour import frame
from f4_formal_cycles import coefficient_matrix, probe
from f4_joined_boundaries import evaluate, native, parent
from mixed33_exact import read


def run(plan, previous, outgoing, work, workers, count, channel="adjoint", source_count=0):
    work.mkdir(exist_ok=True,parents=True)
    z,old,new_out=read(plan),read(previous),read(outgoing)
    prime=old["prime"]
    assert prime==new_out["prime"]
    sources=list(map(m.decode,z["source_dictionary"]))
    before=old["with_formal_image"]
    rows=before["rows"]
    dims={int(d):n for d,n in old["source_character"].items()}
    needed_formal=old["formal_image_word_columns"]
    formal=coefficient_matrix([z["images"][i] for i in rows],len(z["target_dictionary"]),prime)
    src,frames=probe([z["sources"][i] for i in rows],sources,dims,prime,m.library(work,27,prime),
                     workers,507260)
    base=np.column_stack((src,formal[:,needed_formal].toarray()))
    assert int(nmod_mat(base[:,before["columns"]].tolist(),prime).det())==before["determinant"]
    seeds,images=[],[]
    for seed in range(194720,194720+count):
        img=parent(seed,2,channel=channel)
        if img:
            seeds.append(seed)
            images.append(img)
    boundaries=len(images)
    source_seeds=[]
    for seed in range(279331,279331+source_count):
        source=parent(seed,1,channel=channel,source_only=True)
        if source:
            source_seeds.append(seed)
            images.append(source)
    graphs=sorted({g for row in images for g in row})
    grouped={d:[g for g in graphs if m.nferm(g)==d] for d in dims}
    assert sum(map(len,grouped.values()))==len(graphs)
    lib=native(work,prime,channel)
    parts=[]
    for d,gs in sorted(grouped.items()):
        ids={g:i for i,g in enumerate(gs)}
        co=coefficient_matrix([[(ids[g],n) for g,n in row.items() if g in ids] for row in images],
                              len(gs),prime)
        alphabet={l for w in sources if m.nferm(w)==d for tr in w for l in tr}
        for spec in (v for v in frames if v["degree"]==d):
            f=frame(alphabet,spec["points"],spec["k"],d,prime,spec["seed"])
            v=evaluate(gs,f,work,workers,lib,channel)
            parts.append(np.asarray(co@v)%prime)
        print("boundary evaluations",d,len(gs),flush=True)
    values=np.column_stack(parts)
    with_zero=np.column_stack((values[:boundaries],np.zeros((boundaries,len(needed_formal)),np.int64)))
    extended=np.vstack((base,with_zero))
    joint=m.minor(extended,prime)
    lower=joint["rank"]-old["formal_terminal_singlets"]
    kernel_upper=sum(dims.values())-new_out["augmented"]["rank"]
    assert 0 <= lower <= kernel_upper
    record=dict(prime=prime,attempted=count,seeds=seeds,paired_words=len(graphs),
            augmented=joint,added_seeds=[seeds[i-len(rows)] for i in joint["rows"] if i>=len(rows)],
            old_closed_lower=old["restricted_formal_kernel_lower"],closed_and_boundary_lower=lower,
            kernel_upper=kernel_upper,quotient_upper=kernel_upper-lower)
    if source_seeds:
        expanded=m.minor(np.vstack((src,values[boundaries:])),prime)
        assert old["source"]["rank"]<=expanded["rank"]<=sum(dims.values())
        record.update(source_seeds=source_seeds,source_rank_with_new_parents=expanded,
                      paired_boundary_words=len({g for row in images[:boundaries] for g in row}))
    result_path=work/("jordan_incoming.json" if channel=="jordan" else "incoming.json")
    m.write(result_path,record)
    print("joint rank",joint["rank"],"new closed lower",lower,"quotient upper",kernel_upper-lower,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan",type=Path,required=True)
    p.add_argument("--previous",type=Path,required=True)
    p.add_argument("--outgoing",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=16)
    p.add_argument("--count",type=int,default=128)
    p.add_argument("--channel",choices=("adjoint","jordan"),default="adjoint")
    p.add_argument("--source-count",type=int,default=0)
    assert __debug__
    run(**vars(p.parse_args()))
