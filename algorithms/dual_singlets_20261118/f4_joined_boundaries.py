"""Classical singlet boundaries built by pairing two colour covariants in the F4 adjoint."""
import argparse
import ctypes
import math
import random
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from f4_colour import basis, frame, jordan_colours
from f4_formal_cycles import coefficient_matrix, probe
from fortuity.full.letters import Q_letter, is_ferm
from fortuity.hwcoh.evalg import table
from mixed33_exact import read


@lru_cache(None)
def open_word(w, channel):
    n = sum(map(is_ferm,w))
    sign = (-1)**(len(w)+(channel=="adjoint")+n*(n-1)//2)
    rev = w[::-1]
    if rev == w and sign < 0:
        return 0,w
    return (sign,rev) if rev<w else (1,w)


def graph(pair, channel):
    (s,a),(t,b) = (open_word(w,channel) for w in pair)
    n = sum(map(is_ferm,a))*sum(map(is_ferm,b))
    if a==b and n%2:
        return 0,(a,b)
    return (s*t*(-1)**n,(b,a)) if b<a else (s*t,(a,b))


def vary(expr, axis, channel):
    out = Counter()
    for pair,c in expr.items():
        phase = 1
        for j,w in enumerate(pair):
            for k,l in enumerate(w):
                a,b,theta=l
                if axis==0:
                    terms=[(phase*n,(x,y)) for n,x,y in Q_letter(l)]
                elif axis==3:
                    terms=[(a+1,((a+1,b-1,theta),))] if b else []
                else:
                    target=tuple(sorted(axis-1 if x==axis else x for x in theta))
                    terms=[(1,((a,b,target),))] if axis in theta and axis-1 not in theta else []
                for n,letters in terms:
                    new=w[:k]+letters+w[k+1:]
                    sign,g=graph((new,pair[1]) if j==0 else (pair[0],new),channel)
                    if sign:
                        out[g]+=c*n*sign
                if is_ferm(l):
                    phase=-phase
    return {g:c for g,c in out.items() if c}


def parent(seed, deficit=1, channel="adjoint", source_only=False):
    assert channel in ("adjoint","jordan")
    rng=random.Random(seed)
    slots=list(range(11))
    rng.shuffle(slots)
    # Pair or triple derivatives on a superfield before applying its classical variation.
    if deficit==1:
        pattern=[tuple(slots[:2])]+[(x,) for x in slots[2:]]
    elif deficit==2 and seed%2:
        pattern=[tuple(slots[:3])]+[(x,) for x in slots[3:]]
    elif deficit==2:
        pattern=[tuple(slots[:2]),tuple(slots[2:4])]+[(x,) for x in slots[4:]]
    else:
        raise ValueError(deficit)
    rng.shuffle(pattern)
    cut=rng.choice((3,4,5) if channel=="adjoint" else (4,5))
    result=Counter()
    for bits in product((0,1),repeat=4):
        spins=tuple(j for b in bits for j in (b,1-b))
        for flavours in permutations(range(3)):
            n=(-1)**(sum(bits)+sum(i>j for p,i in enumerate(flavours) for j in flavours[p+1:]))
            word=[]
            for spec in pattern:
                dz=[spins[i] for i in spec if i<8]
                th=[flavours[i-8] for i in spec if i>=8]
                a,b=dz.count(0),dz.count(1)
                n*=math.factorial(a)*math.factorial(b)*(-1)**sum(i>j for p,i in enumerate(th) for j in th[p+1:])
                word.append((a,b,tuple(sorted(th))))
            s,g=graph((tuple(word[:cut]),tuple(word[cut:])),channel)
            if s:
                result[g]+=s*n
    source={g:c for g,c in result.items() if c}
    assert all(not vary(source,axis,channel) for axis in (1,2,3))
    if source_only:
        return source
    image=vary(source,0,channel)
    if deficit==2:
        assert not vary(image,0,channel)
    return image


def native(work,prime,channel="adjoint"):
    source=Path(__file__).with_name("f4_colour_contractions.cpp")
    path=work/("open.so" if channel=="adjoint" else "open_jordan.so")
    subprocess.run(["c++","-O3","-std=c++17","-shared","-fPIC","-I",str(m.RT),
                    "-DTRACE_FRAME_MATRIX_SIZE=27",f"-DTRACE_FRAME_PRIME={prime}",
                    f"-DF4_COLOUR_COUNT={len(insertions(channel))}",
                    str(source),"-o",str(path)],check=True)
    lib=ctypes.CDLL(str(path))
    i,p=ctypes.c_int,ctypes.c_void_p
    lib.f4_open.argtypes=[i]*5+[p]*3+[i,p,p,i,p,p]
    lib.f4_open.restype=i
    lib.f4_pair.argtypes=[i]*5+[p]*8
    lib.f4_pair.restype=None
    return lib


def insertions(channel):
    if channel=="jordan":
        return jordan_colours()
    assert channel=="adjoint"
    return basis()[0]


def colour_inputs(f, channel):
    alphabet=sorted(f.L)
    ids={l:i+1 for i,l in enumerate(alphabet)}
    parts=[np.ascontiguousarray(f.L[l][1],dtype=np.int64).ravel() for l in alphabet]
    off=np.asarray([0,0]+list(np.cumsum([len(a) for a in parts])[:-1]),np.uint64)
    parity=np.asarray([0]+[f.L[l][0] for l in alphabet],np.int32)
    mats=insertions(channel)
    gram=np.einsum("iab,jba->ij",mats,mats)
    assert np.array_equal(gram,np.diag(np.diag(gram)))
    terms=np.array([(int(i),int(b*27+a),int(mats[i,a,b])) for i,a,b in np.argwhere(mats)],np.int32)
    inv=np.array([pow(int(x),-1,f.p) for x in np.diag(gram)],np.int64)
    return ids,parity,off,np.concatenate(parts),terms,inv


def open_values(words,degree,f,lib,inputs):
    ids,parity,off,values,gen,_=inputs
    encoded=[bytes(ids[l] for l in w) for w in words]
    starts=np.asarray([0]+list(np.cumsum(list(map(len,encoded)))),np.uint64)
    raw=np.frombuffer(b"".join(encoded),dtype=np.uint8)
    out=np.empty((len(words),f.P,math.comb(f.k,degree),len(inputs[-1])),np.uint32)
    rc=lib.f4_open(f.k,f.P,degree,f.p,len(ids),parity.ctypes.data,off.ctypes.data,values.ctypes.data,
                   len(words),starts.ctypes.data,raw.ctypes.data,len(gen),gen.ctypes.data,out.ctypes.data)
    assert rc==0,rc
    return out


def pair_values(pairs,dl,dr,left,right,f,lib):
    i,j,s,st=table(f.k,dl,dr)
    i,j,s,st=(np.asarray(x,np.int32) for x in (i,j,s,np.r_[st,len(i)]))
    ids=np.asarray(pairs,np.int32)
    out=np.empty((len(pairs),f.P,math.comb(f.k,dl+dr)),np.int64)
    assert left.shape[-1]*int(max(np.diff(st)))*(f.p-1)**2 < 2**63
    lib.f4_pair(len(pairs),f.P,left.shape[2],right.shape[2],out.shape[2],ids.ctypes.data,
                i.ctypes.data,j.ctypes.data,s.ctypes.data,st.ctypes.data,left.ctypes.data,
                right.ctypes.data,out.ctypes.data)
    return out.reshape(len(pairs),-1)


def explicit(w,f,channel):
    d,mat=f.L[w[0]]
    for l in w[1:]:
        mat,d=f.mat_mul_letter(mat,d,l)
    mats=insertions(channel)
    return np.einsum("iab,pxba->pxi",mats,mat)%f.p


def evaluate(graphs,f,work,workers,lib=None,channel="adjoint"):
    if not graphs:
        return np.zeros((0,f.ncols),np.int64)
    if lib is None:
        lib=native(work,f.p,channel)
    inputs=colour_inputs(f,channel)
    words=sorted({w for g in graphs for w in g})
    lists=defaultdict(list)
    for w in words:
        lists[sum(map(is_ferm,w))].append(w)
    records,indices={},{}
    for degree,ws in lists.items():
        indices.update({w:i for i,w in enumerate(ws)})
        chunks=np.array_split(np.arange(len(ws)),min(len(ws),workers))
        with ThreadPoolExecutor(workers) as pool:
            v=np.concatenate(list(pool.map(lambda js:open_values([ws[i] for i in js],degree,f,lib,inputs),chunks)))
        np.testing.assert_array_equal(v[0],explicit(ws[0],f,channel))
        records[degree]=v
    grouped=defaultdict(list)
    for idx,(a,b) in enumerate(graphs):
        grouped[sum(map(is_ferm,a)),sum(map(is_ferm,b))].append((idx,indices[a],indices[b]))
    full=np.empty((len(graphs),f.ncols),np.int64)
    for (dl,dr),items in grouped.items():
        left=(records[dl].astype(np.int64)*inputs[-1]%f.p).astype(np.uint32)
        right=records[dr]
        chunks=np.array_split(np.asarray(items,np.int32), min(len(items),workers))
        def do(chunk):
            v=pair_values(np.ascontiguousarray(chunk[:,1:]),dl,dr,left,right,f,lib)
            full[chunk[:,0]]=v
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(do,chunks))
        idx,a,b=items[0]
        ref=np.zeros((f.P,math.comb(f.k,dl+dr)),np.int64)
        for n in range(len(inputs[-1])):
            val,_=f.scal_mul(left[a,:,:,n].astype(np.int64),dl,right[b,:,:,n].astype(np.int64),dr)
            ref=(ref+val)%f.p
        np.testing.assert_array_equal(full[idx],ref.ravel())
    return full


def run(plan,previous,work,workers,count,channel="adjoint"):
    work.mkdir(exist_ok=True,parents=True)
    z,old=read(plan),read(previous)
    prime=old["prime"]
    chosen=old["classical"]["rows"]
    baseline_rows=[z["images"][i] for i in chosen]
    dct=list(map(m.decode,z["target_dictionary"]))
    needed=sorted({i for row in baseline_rows for i,c in row})
    ix={i:j for j,i in enumerate(needed)}
    base,frames=probe([[(ix[i],v) for i,v in row] for row in baseline_rows],[dct[i] for i in needed],
                      {8:144},prime,m.library(work,27,prime),workers,507100)
    assert len(frames)==1
    minor=base[:,old["classical"]["columns"]]
    assert int(nmod_mat(minor.tolist(),prime).det())==old["classical"]["determinant"]
    seeds,images=[],[]
    for seed in range(279331,279331+count):
        img=parent(seed,channel=channel)
        if img:
            seeds.append(seed)
            images.append(img)
    graphs=sorted({g for row in images for g in row})
    alphabet={l for g in graphs for w in g for l in w}
    fr=frames[0]
    f=frame(alphabet,fr["points"],fr["k"],8,prime,fr["seed"])
    print("graph parents",len(seeds),"graph monomials",len(graphs),"open words",
          len({w for g in graphs for w in g}),flush=True)
    v=evaluate(graphs,f,work,workers,channel=channel)
    ids={g:i for i,g in enumerate(graphs)}
    co=coefficient_matrix([[(ids[g],n) for g,n in row.items()] for row in images],len(graphs),prime)
    new=np.asarray(co@v)%prime
    result=m.minor(np.vstack((base,new)),prime)
    assert old["classical"]["rank"]<=result["rank"]<=144
    result_path=work/("jordan_outgoing.json" if channel=="jordan" else "joined.json")
    m.write(result_path,dict(prime=prime,seeds=seeds,attempted=count,paired_words=len(graphs),
            base_rows=chosen,new_independent_seeds=[seeds[i-len(chosen)] for i in result["rows"] if i>=len(chosen)],
            augmented=result,previous_rank=old["classical"]["rank"],
            kernel_upper=745-result["rank"],quotient_upper=745-result["rank"]-old["restricted_formal_kernel_lower"]))
    print("outgoing rank",result["rank"],"det",result["determinant"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan",type=Path,required=True)
    p.add_argument("--previous",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=16)
    p.add_argument("--count",type=int,default=144)
    p.add_argument("--channel",choices=("adjoint","jordan"),default="adjoint")
    assert __debug__
    run(**vars(p.parse_args()))
