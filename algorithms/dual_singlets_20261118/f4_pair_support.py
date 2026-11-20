"""Project the classical closed span to the only fermion degree seen by Q1."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np

import mixed_incoming as m
from f4_colour import frame
from f4_formal_cycles import coefficient_matrix, probe
from f4_joined_boundaries import evaluate, native, parent
from mixed33_exact import read
from rank_four_s5 import Q1d_multi
from q1d_exact import pair_terms


def run(plan, previous, jordan, out, workers):
    z,old,j=map(read,(plan,previous,jordan))
    p=old["prime"]
    assert j["prime"]==p
    out.mkdir(exist_ok=True,parents=True)
    words=list(map(m.decode,z["source_dictionary"]))
    tally=Counter()
    for w in words:
        letters=[l for name in w for l in name]
        flavours=[l for l in letters if l[2]]
        pairs=[(a,b) for k,a in enumerate(flavours) for b in flavours[k+1:] if pair_terms(a,b)]
        degree=m.nferm(w)
        assert (degree,len(flavours),len(pairs)) in ((7,3,0),(9,2,1))
        if pairs:
            assert all(not l[0] and not l[1] for l in flavours)
        tally[degree]+=1
    # Word indices, including words that never appear in a chosen protected row.
    chosen=[i for i,w in enumerate(words) if m.nferm(w)==9]
    remap={i:k for k,i in enumerate(chosen)}
    rows=[[(remap[i],c) for i,c in row if i in remap] for row in z["sources"]]
    ws=[words[i] for i in chosen]
    src,frames=probe(rows,ws,{9:79},p,m.library(out,27,p),workers,507260)
    formal=coefficient_matrix(z["images"],len(z["target_dictionary"]),p)[:,old["formal_image_word_columns"]].toarray()
    one=m.minor(src,p)
    with_formal=np.column_stack((src,formal))
    old_joint=m.minor(with_formal,p)
    boundaries=[parent(i,2,"jordan") for i in j["added_seeds"]]
    raw=[parent(i,1,"jordan",source_only=True) for i in j["source_seeds"]]
    recipes=boundaries+raw
    graphs=sorted({g for row in recipes for g in row if m.nferm(g)==9})
    ids={g:i for i,g in enumerate(graphs)}
    co=coefficient_matrix([[(ids[g],c) for g,c in row.items() if g in ids] for row in recipes],len(graphs),p)
    alphabet={l for w in ws for tr in w for l in tr}
    assert all(l in alphabet for g in graphs for w in g for l in w)
    vals=[]
    lib=native(out,p,"jordan")
    for spec in frames:
        f=frame(alphabet,spec["points"],spec["k"],9,p,spec["seed"])
        vals.append(np.asarray(co@evaluate(graphs,f,out,workers,lib,"jordan"))%p)
    coeff=np.column_stack(vals)
    b=len(boundaries)
    extended=np.vstack((with_formal,np.column_stack((coeff[:b],np.zeros((b,formal.shape[1]),np.int64)))))
    joint=m.minor(extended,p)
    proj=joint["rank"]-old["formal_terminal_singlets"]
    full=m.minor(np.vstack((src,coeff)),p)
    assert 0<=proj<=full["rank"]<=79
    record=dict(prime=p,weight=z["mu"],word_counts=dict(tally),ambient_rank=one,augmented_ambient_rank=full,
        formal_rank=old["formal_terminal_singlets"],formal_joint=old_joint,joint_with_boundaries=joint,
        protected_projection_lower=proj,unaccounted_projected_upper=79-proj,
        added_boundary_seeds=j["added_seeds"],raw_source_seeds=j["source_seeds"],frames=frames)
    m.write(out/"pair_support.json",record)
    print("nine-fermion projection",one["rank"],full["rank"],"closed lower",proj,"remaining upper",79-proj,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("plan","previous","jordan","out"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--workers",type=int,default=8)
    assert __debug__
    run(**vars(p.parse_args()))
