"""Compute filtered maps for the SO(7)/USp(6) level-24 non-Coulomb comparison."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from flint import nmod_mat
from scipy.sparse import csr_matrix

import mixed_incoming as m
from filtered_rank import filtered_ranks
from trace_descent import action, dimensions, highest_weights, settings


def dense_frame(colours, group, points, degree, alphabet, seed):
    return m.Frame(colours//2, points, degree+1, [degree], alphabet, seed=seed, p=m.P,
                   group="SO_odd" if group=="SO" else "Sp")


def coordinates(rows, degree, alphabet, samples, colours, group, library, workers):
    words=sorted({w for row in rows for w in row if m.nferm(w)==degree})
    width=sum(len(columns) for _,_,columns in samples)
    if not words:
        return np.zeros((len(rows),width),dtype=np.int64)
    ids={w:i for i,w in enumerate(words)}
    entries,indices,starts=[],[],[0]
    for row in rows:
        for w,n in row.items():
            if w in ids:
                n=Fraction(n)
                entries.append(n.numerator*pow(n.denominator,-1,m.P)%m.P)
                indices.append(ids[w])
        starts.append(len(entries))
    assert max(np.diff(starts),default=0)*(m.P-1)**2 < 2**63
    coefficients=csr_matrix((np.asarray(entries,np.int64),indices,starts),shape=(len(rows),len(words)))
    def evaluate_sample(sample):
        points,seed,keep=sample
        frame=dense_frame(colours,group,points,degree,alphabet,seed)
        result=np.zeros((len(rows),len(keep)),dtype=np.int64)
        for lo in range(0,len(words),3072):
            hi=min(lo+3072,len(words))
            values=m.values(words[lo:hi],frame,library)
            result=(result+coefficients[:,lo:hi]@values[:,keep])%m.P
        return result
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return np.column_stack(list(pool.map(evaluate_sample,samples)))


def character_frame(rows,dim,degree,alphabet,colours,group,library,workers,seed):
    points=math.ceil((dim+12)/(degree+1))
    samples=[(min(8,points-i),seed+6151*i,list(range(min(8,points-i)*(degree+1))))
             for i in range(0,points,8)]
    table=coordinates(rows,degree,alphabet,samples,colours,group,library,workers)
    for k in range(8):
        field=nmod_mat(table.tolist(),m.P)
        columns=m.pivotcols(field)
        if len(columns)==dim:
            break
        assert len(columns)<dim
        extra=(8,seed+6151*(points+8*k),list(range(8*(degree+1))))
        samples.append(extra)
        table=np.column_stack((table,coordinates(rows,degree,alphabet,[extra],colours,group,library,workers)))
    else:
        raise ValueError(("incomplete target frame",group,degree,dim,len(columns)))
    piv=m.pivotcols(nmod_mat(table[:,columns].T.tolist(),m.P))
    det=int(nmod_mat(table[np.ix_(piv,columns)].tolist(),m.P).det())
    assert det and len(piv)==dim
    kept=set(columns)
    compact=[]
    offset=0
    for n,point_seed,jj in samples:
        selected=[j for j in jj if offset+j in kept]
        if selected: compact.append((n,point_seed,selected))
        offset+=len(jj)
    return compact,dict(character_dimension=dim,formal_rows=len(rows),minor_rows=piv,
                         determinant_mod_prime=det,samples=compact)


def weight(grade,flavours):
    return (grade,grade)+tuple(n-grade for n in flavours)


def run_grade(flavours,deficit,grade,group,colours,out,workers):
    settings()
    out.mkdir(parents=True,exist_ok=True)
    w=weight(grade,flavours)
    end=sum(w)-deficit
    target=highest_weights(w,end)
    dt=dimensions(w,end,colours,group)
    alphabet=sorted({l for content in m.contents(w,end) for l in content})
    library=m.library(out,colours,m.P)
    info=dict(flavours=flavours,group=group,colours=colours,grade=grade,deficit=deficit,
              weight=w,target_length=end,target_dimensions=dt,prime=m.P,frames={},blocks={})
    frames={}
    for degree,dim in dt.items():
        frames[degree],info["frames"][degree]=character_frame(target[degree],dim,degree,alphabet,
            colours,group,library,workers,192107+971*grade+13*degree+colours)
        print(group,deficit,grade,"frame",degree,dim,flush=True)
    maximum=(deficit+2)//2
    for source_grade in reversed(range(max(0,grade-maximum),grade+1)):
        source_weight=weight(source_grade,flavours)
        length=sum(source_weight)-deficit-1
        src=highest_weights(source_weight,length)
        source=[row for deg in sorted(src) for row in src[deg]]
        order=grade-source_grade
        cells=np.empty((len(source),sum(dt.values())),dtype=np.int64)
        for i in range(0,len(source),384):
            actions=[action(row,order,colours,group) for row in source[i:i+384]]
            if dt:
                cells[i:i+len(actions)]=np.column_stack([
                    coordinates(actions,degree,alphabet,frame,colours,group,library,workers)
                    for degree,frame in frames.items()])
            else: assert not any(actions) or not sum(dt.values())
        filename=f"block{source_grade}_{grade}.npy"
        stage=out/(filename+".part")
        with stage.open("wb") as handle: np.save(handle,cells,allow_pickle=False)
        stage.replace(out/filename)
        info["blocks"][filename]=dict(shape=cells.shape,sha256=hashlib.sha256((out/filename).read_bytes()).hexdigest())
        m.write(out/"progress.json",info)
        print(group,deficit,grade,"Q",order,cells.shape,flush=True)
    m.write(out/"grade.json",info)
    return info


def assemble(flavours,deficit,group,colours,out):
    source=[];target=[];blocks={};hashes={}
    for grade in range(min(flavours)+1):
        folder=out/f"g{grade}"
        record=json.loads((folder/"grade.json").read_bytes())
        assert (record["flavours"],record["group"],record["colours"],record["deficit"],record["grade"])==(
            list(flavours),group,colours,deficit,grade)
        for name,values in record["blocks"].items():
            a,b=map(int,name.removeprefix("block").removesuffix(".npy").split("_"))
            assert b==grade and hashlib.sha256((folder/name).read_bytes()).hexdigest()==values["sha256"]
            blocks[a,b]=np.load(folder/name,allow_pickle=False)
        s,t=blocks[grade,grade].shape;source.append(s);target.append(t)
        hashes[str(grade)]=hashlib.sha256((folder/"grade.json").read_bytes()).hexdigest()
    assert all(v.shape==(source[a],target[b]) for (a,b),v in blocks.items())
    intervals=filtered_ranks(blocks,source,target,m.P)
    result=dict(flavours=flavours,group=group,colours=colours,target_deficit=deficit,prime=m.P,
                formal_source_dimensions=source,target_dimensions=target,grade_sha256=hashes,
                interval_ranks_and_pairings=intervals,
                ranks_by_order=[sum(n for a,b,r,n in intervals if b-a==k) for k in range(len(source))])
    m.write(out/"band.json",result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flavours",type=int,nargs=3,default=(6,3,3))
    parser.add_argument("--group",choices=("SO","Sp"),required=True)
    parser.add_argument("--colours",type=int,required=True)
    parser.add_argument("--deficit",type=int,required=True,choices=(0,1,2))
    parser.add_argument("--grade",type=int)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    assert __debug__ and ((args.group=="SO" and args.colours%2) or (args.group=="Sp" and args.colours%2==0))
    if args.grade is None: print(assemble(args.flavours,args.deficit,args.group,args.colours,args.out)["ranks_by_order"])
    else: run_grade(args.flavours,args.deficit,args.grade,args.group,args.colours,args.out/f"g{args.grade}",args.workers)
