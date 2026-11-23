"""Exact singlet sources and finite-colour incoming Q0 minors at bidegree (3,3).

The source kernel is computed over Z, separately in invariant trace-shape blocks.
Matrices of trace evaluations are rectangular witnesses, never a priori bases.
"""
import argparse, ctypes, gzip, hashlib, json, math, os, subprocess, sys, tempfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from flint import fmpz_mat, nmod_mat
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
DEFAULT_RT = HERE.parent/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
RT = Path(os.environ.get('N4_TRACE_RUNTIME', DEFAULT_RT)).resolve()
sys.path.insert(0, str(RT))
import sp_environment
from g2_full_classical_screen_v1 import contents, sector
from full_hw_so import so_canon, shape, raise_word
from fortuity.full import words as fw
from fortuity.full.evalg import Frame, nferm
from fortuity.full.wordenum import multitraces_with_content

MU=(3,3,1,1,1)
P=1000003

def digest(b):return hashlib.sha256(b).hexdigest()
def compact(v):return json.dumps(v,separators=(',',':')).encode()
def encode(w):return [[[i,j,list(t)] for i,j,t in tr] for tr in w]
def write(path,obj):
    data=compact(obj)+b'\n'
    if path.suffix=='.gz': data=gzip.compress(data,mtime=0)
    temp=path.with_suffix(path.suffix+'.part');temp.write_bytes(data);temp.replace(path)

def words_by_shape(length,mu):
    blocks=defaultdict(set)
    for content in contents(mu,length):
        for w in multitraces_with_content(content):
            s,v=so_canon(w)
            if s:blocks[shape(v,'spin')].add(v)
    return {k:sorted(v) for k,v in sorted(blocks.items())}

@lru_cache(maxsize=65536)
def differential(w):
    out=Counter()
    for tr,c in fw.Q_multi(w).items():
        s,v=so_canon(tr)
        if s:out[v]+=s*c
    return tuple((v,c) for v,c in out.items() if c)

def q(expr):
    out=Counter()
    for w,n in expr.items():
        for v,c in differential(w):out[v]+=n*c
    return {v:c for v,c in out.items() if c}

def singlets(words):
    # All simple raising operators preserve this cyclic/reversal trace shape.
    equations=[{(i,v):n for i in (1,2,3) for v,n in raise_word(w,i).items()} for w in words]
    ids=sorted(set().union(*(e.keys() for e in equations)));rows={key:i for i,key in enumerate(ids)}
    a=fmpz_mat(len(ids),len(words))
    for j,e in enumerate(equations):
        for k,c in e.items():a[rows[k],j]=c
    z,dim=a.nullspace()
    for j in range(dim):
        co=[int(z[i,j]) for i in range(len(words))];g=math.gcd(*co)
        if next(c for c in co if c)<0:g=-g
        co=[c//g for c in co]; assert (a*fmpz_mat([[c] for c in co])).is_zero()
        yield {w:c for w,c in zip(words,co) if c}

def flatten(rows):
    words=sorted(set().union(*(e.keys() for e in rows)));wid={w:i for i,w in enumerate(words)}
    terms=[sorted((wid[w],c) for w,c in e.items()) for e in rows];return words,terms

def prepare(mu,length):
    blocks=words_by_shape(length,mu); sources=[];images=[]
    for block in blocks.values():
        for e in singlets(block):
            d=q(e);assert not q(d)
            sources.append(e);images.append(d)
    words,rows=flatten(images);srcwords,sc=flatten(sources)
    return dict(mu=mu,length=length,source_words=sum(map(len,blocks.values())),source_shapes=len(blocks),formal_singlets=len(rows),source_sector=sector(mu,length),target_sector=sector(mu,length+1),source_dictionary=[encode(w) for w in srcwords],sources=sc,target_dictionary=[encode(w) for w in words],images=rows)

def decode(w):return tuple(tuple((i,j,tuple(t)) for i,j,t in tr) for tr in w)

def library(dirpath,n,p):
    src=RT/'trace_frame_bmn_rank_threadsafe_v1.cpp';target=dirpath/f'frame_{n}_{p}.so'
    subprocess.run(['c++','-O3','-std=c++17','-shared','-fPIC',f'-DTRACE_FRAME_MATRIX_SIZE={n}',f'-DTRACE_FRAME_PRIME={p}',str(src),'-o',str(target)],check=True)
    lib=ctypes.CDLL(str(target));lib.trace_frame_error.restype=ctypes.c_char_p
    lib.trace_frame_evaluate.argtypes=[ctypes.c_int]*5+[ctypes.c_void_p]*3+[ctypes.c_int]+[ctypes.c_void_p]*3
    lib.trace_frame_evaluate.restype=ctypes.c_int;return lib

def values(words,frame,lib):
    ids={l:i for i,l in enumerate(sorted(frame.L),1)};pars=[0];offs=[0];arr=[];n=0
    for l in sorted(ids):
        d,v=frame.L[l];pars.append(d);offs.append(n);a=np.ascontiguousarray(v,dtype=np.int64).reshape(-1);n+=len(a);arr.append(a)
    pars=np.asarray(pars,np.int32);offs=np.asarray(offs,np.uint64);mat=np.concatenate(arr);raw=bytearray(); starts=[0]
    for w in words:
        for t in w:raw.extend(ids[l] for l in t); raw.append(0)
        starts.append(len(raw))
    starts=np.asarray(starts,np.uint64);raw=np.frombuffer(bytes(raw),np.uint8);out=np.empty((len(words),frame.ncols),np.int64)
    ok=lib.trace_frame_evaluate(frame.k,frame.P,frame.degrees[0],frame.p,len(ids),pars.ctypes.data,offs.ctypes.data,mat.ctypes.data,len(words),starts.ctypes.data,raw.ctypes.data,out.ctypes.data)
    assert ok==0,lib.trace_frame_error();return out

def pivotcols(a):
    r,k=a.rref();p=[];j=0
    for i in range(k):
        while not r[i,j]: j+=1
        p.append(j);j+=1
    return p

def minor(values,p):
    a=nmod_mat(values.tolist(),p);cols=pivotcols(a);rows=pivotcols(a.transpose());assert len(rows)==len(cols)
    sq=values[np.ix_(rows,cols)];det=int(nmod_mat(sq.tolist(),p).det());assert det
    return dict(rank=len(cols),rows=rows,columns=cols,determinant=det,minor_sha256=digest(sq.astype('<u4').tobytes()))

def run(plan,group,prime,work,workers,pts,k,seed,check_count):
    words=list(map(decode,plan['target_dictionary']));rows=plan['images']; degrees=sorted(set(map(nferm,words))); take=defaultdict(list)
    for i,w in enumerate(words):take[nferm(w)].append(i)
    assert max(sum(c%prime for _,c in e) for e in rows)*(prime-1)<2**63
    rr,cc,coeff=zip(*((i,j,c%prime) for i,e in enumerate(rows) for j,c in e));a=csr_matrix((np.asarray(coeff,np.int64),(rr,cc)),shape=(len(rows),len(words)))
    lib=library(work,7 if group=='SO_odd' else 6,prime);by_degree=[];fs=[]
    assert len(pts)==len(degrees) and max(degrees)<=k
    for d,points in zip(degrees,pts):
        indices=take[d];ws=[words[i] for i in indices];alphabet={l for w in ws for t in w for l in t};f=Frame(3,points,k,[d],alphabet,seed=seed+d,p=prime,group=group)
        chunks=np.array_split(np.arange(len(ws)),min(len(ws),workers));
        with ThreadPoolExecutor(max_workers=workers) as ex:
            bs=list(ex.map(lambda sub:values([ws[i] for i in sub],f,lib),chunks))
        ev=np.concatenate(bs);del bs
        control=sorted(set(np.linspace(0,len(ws)-1,min(check_count,len(ws)),dtype=int).tolist()));
        if control:np.testing.assert_array_equal(ev[control],f.eval_words([ws[j] for j in control],dtype=np.int64))
        by_degree.append(np.asarray(a[:,indices]@ev)%prime);fs.append(dict(fermions=d,k=k,points=points,seed=seed+d,words=len(ws),columns=f.ncols,independent_python_words=len(control)));del ev
    v=np.concatenate(by_degree,axis=1);m=minor(v,prime)
    return dict(group=group,prime=prime,source_words=plan['source_words'],formal_singlets=len(rows),target_words=len(words),matrix_shape=v.shape,frames=fs,evaluation_sha256=digest(v.astype('<u4').tobytes()),**m)

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--plan',type=Path,required=True);a.add_argument('--prepare',action='store_true');a.add_argument('--mu',default='3,3,1,1,1');a.add_argument('--length',type=int,default=7);a.add_argument('--group',choices=['SO_odd','Sp']);a.add_argument('--output',type=Path);a.add_argument('--prime',type=int,default=P);a.add_argument('--workers',type=int,default=4);a.add_argument('--points',default='6,8');a.add_argument('--k',type=int,default=8);a.add_argument('--seed',type=int,default=3711902);a.add_argument('--check-words',type=int,default=8);args=a.parse_args();assert __debug__
    if args.prepare:
        z=prepare(tuple(map(int,args.mu.split(','))),args.length);write(args.plan,z);print('source words, singlets, target words',z['source_words'],z['formal_singlets'],len(z['target_dictionary']),flush=True)
    if args.group:
        plan=gzip.decompress(args.plan.read_bytes()); z=json.loads(plan)
        with tempfile.TemporaryDirectory() as td:r=run(z,args.group,args.prime,Path(td),args.workers,list(map(int,args.points.split(','))),args.k,args.seed,args.check_words)
        r.update(plan_sha256=digest(plan),evaluation_source_sha256=digest((RT/'trace_frame_bmn_rank_threadsafe_v1.cpp').read_bytes()))
        write(args.output,r);print(r['group'],r['matrix_shape'],r['rank'],r['determinant'],flush=True)
