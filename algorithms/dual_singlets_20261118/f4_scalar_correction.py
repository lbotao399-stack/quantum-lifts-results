"""Test the classical seven-fermion equation for the scalar--matter fifth-power source."""
import argparse
import copy
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

import mixed_incoming as m
from f4_colour import frame
from f4_first_source import singlet, multiply
from f4_formal_cycles import coefficient_matrix, probe
from f4_joined_boundaries import evaluate, native, parent, vary
from fortuity.hwcoh.evalg import table
from mixed33_exact import read
from mixed_middle import pin


def spectator_values(graphs, fr, work, workers, lib, channel="jordan"):
    out=np.empty((len(graphs),fr.ncols),np.int64)
    traces=fr.necklace_values([g[2] for g in graphs])
    groups=defaultdict(list)
    for i,(a,b,tr) in enumerate(graphs):
        groups[m.nferm((a,b))].append((i,(a,b),tr))
    for degree,items in groups.items():
        sf=copy.copy(fr)
        cpp=math.comb(sf.k,degree)
        sf.ncols,sf.cpp,sf.degrees,sf.ncomp,sf.offsets=sf.P*cpp,cpp,[degree],{degree:cpp},{degree:0}
        pairs=sorted({pair for _,pair,_ in items})
        ids={pair:i for i,pair in enumerate(pairs)}
        pair_rows=evaluate(pairs,sf,work,workers,lib,channel).reshape(len(pairs),sf.P,cpp)
        dtrace=fr.degrees[0]-degree
        left,right,sign,starts=table(sf.k,degree,dtrace)
        assert max(np.diff(np.r_[starts,len(left)]))*(sf.p-1)**2<2**63
        for start in range(0,len(items),48):
            chunk=items[start:start+48]
            a=pair_rows[[ids[g] for _,g,_ in chunk]]
            b=np.stack([traces[t][0] for _,_,t in chunk])
            value=np.add.reduceat(a[:,:,left]*b[:,:,right]*sign,starts,axis=2)%sf.p
            out[[i for i,_,_ in chunk]]=value.reshape(len(chunk),-1)
    return out


def evaluate_rows(recipes, alphabet, frames, p, work, workers, lib, spectator, channel="jordan"):
    graphs=sorted({g for row in recipes for g in row})
    ids={g:i for i,g in enumerate(graphs)}
    co=coefficient_matrix([[(ids[g],v) for g,v in row.items()] for row in recipes],len(graphs),p)
    columns=[]
    for s in frames:
        fr=frame(alphabet,s["points"],s["k"],s["degree"],p,s["seed"])
        if spectator:
            evaluated=spectator_values(graphs,fr,work,workers,lib,channel)
        else:
            evaluated=evaluate(graphs,fr,work,workers,lib,channel)
        columns.append(np.asarray(co@evaluated)%p)
    return np.column_stack(columns),len(graphs)


def run(plan, previous, work, workers, count):
    work.mkdir(exist_ok=True,parents=True)
    z,old=read(plan),read(previous)
    p=old["prime"]
    assert old["plan_sha256"]==pin(plan)
    sword=list(map(m.decode,z["source_dictionary"]))
    tword=list(map(m.decode,z["target_dictionary"]))
    idx=[i for i,row in enumerate(z["sources"]) if m.nferm(sword[row[0][0]])==7]
    assert all(m.nferm(sword[j])==7 for i in idx for j,_ in z["sources"][i])
    used=[i for i,w in enumerate(sword) if m.nferm(w)==7]
    remap={j:i for i,j in enumerate(used)}
    srcwords=[sword[i] for i in used]
    library=m.library(work,27,p)
    q=m.q(multiply(singlet(),4))
    ti={w:j for j,w in enumerate(tword)}
    base_out,target_frames=probe([z["images"][i] for i in idx]+[[(ti[w],n) for w,n in q.items()]],
                                 tword,{8:144},p,library,workers,507100)
    base_source,source_frames=probe([[(remap[j],n) for j,n in z["sources"][i]] for i in idx],
                                 srcwords,{7:666},p,library,workers,507260)
    rhs=base_out[-1]
    images,src=base_out[:-1],base_source
    target_letters={l for w in tword for tr in w for l in tr}
    source_letters={l for w in srcwords for tr in w for l in tr}
    lib=native(work,p,"jordan")
    stages=[]
    labels=list(idx)
    def ranks(name):
        source_rank=m.minor(src,p)
        image=m.minor(images,p)
        augmented=m.minor(np.vstack((images,rhs)),p)
        print(name,"source",source_rank["rank"],"image",image["rank"],"with RHS",augmented["rank"],flush=True)
        return dict(source=source_rank,image=image,augmented=augmented)
    stages.append(dict(family="ordinary trace",ranks=ranks("traces")))
    for spectator,start in ((False,279331),(True,360180)):
        seeds,rows,actions=[],[],[]
        for seed in range(start,start+count):
            source=parent(seed,1,"jordan",source_only=True,spectator=spectator)
            if source and m.nferm(next(iter(source)))==7:
                assert {m.nferm(g) for g in source}=={7}
                seeds.append(seed)
                rows.append(source)
                actions.append(vary(source,0,"jordan"))
        new,words=evaluate_rows(rows,source_letters,source_frames,p,work,workers,lib,spectator)
        dq,qwords=evaluate_rows(actions,target_letters,target_frames,p,work,workers,lib,spectator)
        images=np.vstack((images,dq))
        src=np.vstack((src,new))
        labels.extend(seeds)
        name="Jordan pair times a trace" if spectator else "Jordan pair"
        stages.append(dict(family=name,seed_start=start,attempted=count,seeds=seeds,
                           source_words=words,image_words=qwords,ranks=ranks(name)))
        m.write(work/"scalar_correction.json",dict(prime=p,weight=z["mu"],ordinary_rows=idx,
                formal_plan_sha256=pin(plan),source_frames=source_frames,target_frames=target_frames,stages=stages))
    assert src.shape[0]==images.shape[0]==len(labels)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("plan","previous","work"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--count",type=int,default=400)
    p.add_argument("--workers",type=int,default=16)
    assert __debug__
    run(**vars(p.parse_args()))
