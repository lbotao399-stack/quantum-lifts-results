"""Test the first correction of the level-26 SO(7) product in finite colour."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from dual633_pairing import integral, load_source, target_frame
from dual_charge_bands import coordinates
from f4_bmn_primary import translated
from fortuity.bmn.words import w_gen
from gaugino_primary_order import STEPS, charge
from so8_level24_descendant import raising_images, substitute
from trace_descent import action, dimensions, evaluate, highest_weights, settings


def product(left, right):
    out = Counter()
    for a,n in left.items():
        for b,k in right.items():
            sign,w = m.so_canon(a+b)
            if sign:
                out[w] += sign*n*k
    return {w:n for w,n in out.items() if n}


def lower(row):
    for spin,flavour in STEPS:
        row = substitute(row,spin,flavour,False)
    return row


def source_data(path):
    f = load_source(path)
    parent_image, _ = integral(m.q(f))
    frame = target_frame((3,3,1,1,1),9,771253)
    assert not any(evaluate([parent_image],frame,7,"SO")[0])
    w = translated(w_gen((0,),1))
    assert not m.q(w)
    source = product(lower(f),w)
    assert charge(source) == ((0,0,5,4,4),10)
    assert all(not any(r.values()) for r in raising_images(source))
    closure = m.q(source)
    assert closure == product(lower(m.q(f)),w)
    first = action(source,1,7,"SO")
    assert first and charge(first)==((1,1,4,3,3),10)
    return source,closure,first


def grade_probes(dimension,degree,seed):
    count = math.ceil((dimension+12)/(degree+1))
    return [(min(6,count-i),seed+6151*i,list(range(min(6,count-i)*(degree+1))))
            for i in range(0,count,6)]


def run(path,out,workers):
    settings()
    out.mkdir(exist_ok=True,parents=True)
    w = (1,1,4,3,3)
    size = 10
    source, closure, first = source_data(path)
    nums,den = integral(first)
    target_dims = dimensions(w,size,7,"SO")
    record=dict(weight=w,level=26,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                source_terms=len(source),propagated_classical_terms=len(closure),
                first_terms=len(first),first_denominator=den,target_dimensions=target_dims,
                first_fermion_degrees=sorted({m.nferm(tr) for tr in first}),prime=m.P)
    m.write(out/"source.json",record)
    print("source",record,flush=True)
    parents = highest_weights(w,size-1)
    formal = [row for degree in sorted(parents) for row in parents[degree]]
    letters = sorted({l for content in m.contents(w,size) for l in content})
    libraries=m.library(out,7,m.P)
    # A rank increase on any set of probes separates the image over this field.
    # A matching rank is recorded with the probes for a fresh holdout.
    blocks=[]
    for degree,dim in target_dims.items():
        samples=grade_probes(dim,degree,168227+degree*1223)
        target_file=out/f"classical_{degree}.npy"
        rhs=coordinates([nums],degree,letters,samples,7,"SO",libraries,workers)[0]
        cols=len(rhs)
        if not target_file.exists():
            temp=target_file.with_suffix(".part")
            values=np.lib.format.open_memmap(temp,mode="w+",dtype=np.int64,shape=(len(formal)+1,cols))
            values[-1]=rhs*pow(den,-1,m.P)%m.P
            for start in range(0,len(formal),256):
                images=[m.q(r) for r in formal[start:start+256]]
                values[start:start+len(images)]=coordinates(images,degree,letters,samples,7,"SO",libraries,workers)
                if start%4096==0: print("boundary rows",degree,start,flush=True)
            values.flush()
            del values
            temp.replace(target_file)
        block=np.load(target_file,mmap_mode="r",allow_pickle=False)
        assert block.shape==(len(formal)+1,cols) and np.array_equal(block[-1],rhs*pow(den,-1,m.P)%m.P)
        blocks.append(block)
        m.write(out/f"grade_{degree}.json",dict(degree=degree,character_dimension=dim,
                rows=len(formal),columns=cols,seed=samples[0][1],
                rhs_nonzero=int(np.count_nonzero(rhs))))
    # Build the transpose in small row chunks to avoid a matrix-sized Python list.
    equations=nmod_mat(sum(b.shape[1] for b in blocks),len(formal)+1,m.P)
    offset=0
    for block in blocks:
        for j,row in enumerate(block.T):
            for i,co in enumerate(row):
                equations[offset+j,i]=int(co)
        offset+=block.shape[1]
    reduced,rank=equations.rref(inplace=True)
    pivots=[]
    j=0
    for i in range(rank):
        while not reduced[i,j]: j+=1
        pivots.append(j)
        j+=1
    consistent=not pivots or pivots[-1]<len(formal)
    result=dict(**record,formal_rows=len(formal),probe_count=offset,
                augmented_rank=rank,classical_rank=rank-int(not consistent),
                first_correction_consistent_on_probes=consistent)
    if consistent:
        coeff=[0]*len(formal)
        for i,j in enumerate(pivots): coeff[j]=int(reduced[i,len(formal)])
        v=np.array(coeff,dtype=np.int64)
        assert len(formal)*(m.P-1)**2<2**63
        assert all(np.array_equal(v@b[:-1]%m.P,b[-1]) for b in blocks)
        result["boundary_expression_mod_prime"]=[[i,n] for i,n in enumerate(coeff) if n]
        hold=[]
        image=Counter()
        for i,n in result["boundary_expression_mod_prime"]:
            for word,c in m.q(formal[i]).items(): image[word]+=c*n
        residual={word:(int(n)*pow(den,-1,m.P)-image.pop(word,0))%m.P for word,n in nums.items()}
        residual.update((word,-n%m.P) for word,n in image.items())
        residual={word:n for word,n in residual.items() if n}
        for degree in target_dims:
            samples=[(5,958391+821*degree,list(range(5*(degree+1))))]
            values=coordinates([residual],degree,letters,samples,7,"SO",libraries,workers)[0]
            hold.extend(map(int,values))
        result["holdout_coordinates"]=len(hold)
        result["holdout_nonzero"]=sum(bool(v) for v in hold)
    m.write(out/"first.json",result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--workers",type=int,default=6)
    args=parser.parse_args()
    assert __debug__
    run(args.source,args.out,args.workers)
