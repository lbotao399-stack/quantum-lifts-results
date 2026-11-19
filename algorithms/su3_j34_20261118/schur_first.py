"""Adapt the source basis to the onto projection of Q0 onto first targets 1 and 5.

Matrices have recipes as rows.  A correction is extended using its original
full trace recipe.  No component of Q0 is discarded in the source change.
"""
import argparse, gzip, hashlib, json, math, os, random, sys, time
from concurrent.futures import ProcessPoolExecutor
from itertools import chain
from pathlib import Path
import numpy as np
from flint import nmod_mat
import projected_q0 as maps

P=1000003
INPUT='fd22fbcd0ae0ada903acca1831b3ad70dc17fbb9570c833a8dd96544d01040e1'
PLAN='1b7ac3dec77fcb50473d01176a8f309773e8411f0fd4e969f58299a026df6e8c'
SOURCES=(1,2,4,10)
TARGETS=(1,5)
OTHERS=(0,2,3,4,6,7,8,9,11,12)

def digest(x):return hashlib.sha256(x).hexdigest()
def fsha(x):
    with Path(x).open('rb') as f:
        h=hashlib.sha256()
        for part in iter(lambda:f.read(1048576),b''):h.update(part)
        return h.hexdigest()
def pack(a):return np.ascontiguousarray(a,dtype='<u4')
def ahash(a):return digest(pack(a).tobytes())
def save_json(path,data):
    b=(json.dumps(data,indent=2)+'\n').encode();tmp=path.with_suffix('.tmp');tmp.write_bytes(b);tmp.replace(path)
def save_z(path,data):
    b=gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0);tmp=path.with_suffix('.tmp');tmp.write_bytes(b);tmp.replace(path)
def read_inputs(planfile,minor):
    assert fsha(planfile)==PLAN and fsha(minor)==INPUT
    p=json.loads(planfile.read_bytes())['incoming'];c=json.loads(gzip.decompress(minor.read_bytes()));r=c['map']
    assert r==dict(order=1,prime=P,sources=list(SOURCES),targets=list(TARGETS),plan_sha256=PLAN,frame_width=1131,source_dimension=5012,target_dimension=1112)
    assert set(TARGETS).isdisjoint(map(int,p['first_rhs_trace_terms']))
    assert SOURCES==tuple(s['index'] for s in p['sources'] if set(s['targets']).intersection(TARGETS))
    assert OTHERS==tuple(sorted({t for s in SOURCES for t in p['sources'][s]['targets']}-set(TARGETS)))
    words=[(i,maps.decode(w)) for i,w in c['rows']];cols=c['columns'];assert len(words)==len(set(words))==len(cols)==len(set(cols))==1112 and c['determinant']==46239
    return p,words,cols

def batches(words,size):return (words[i:i+size] for i in range(0,len(words),size))

def make_factor(p,words,cols,a):
    from exact_blas_rows import matmul_mod
    maps.setup(a.runtime.resolve(),p,list(TARGETS),list(SOURCES),P);allrows=np.vstack([maps.values(x) for x in batches(words,8)]);A=allrows[:,cols];assert A.shape==(1112,1112)
    m=nmod_mat(A.tolist(),P);det=int(m.det());assert det==46239
    inv=pack(m.inv().tolist());assert np.array_equal(matmul_mod(inv,A,P),np.eye(len(cols),dtype=np.int64));assert np.array_equal(matmul_mod(A,inv,P),np.eye(len(cols),dtype=np.int64))
    fp=a.out/'target_1_5_inverse.npy';np.save(fp,inv,allow_pickle=False)
    r=dict(prime=P,targets=list(TARGETS),sizes=[len(words),allrows.shape[1]],determinant=det,minor_sha256=INPUT,plan_sha256=PLAN,minor_entries_sha256=ahash(A),full_projected_entries_sha256=ahash(allrows),inverse_file_sha256=fsha(fp),both_inverse_products_identity=True,other_targets=list(OTHERS),zero_rhs=True);save_json(a.out/'factor.json',r);print(json.dumps(r),flush=True)

def basis_frame(p,src):
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents
    s=p['sources'][src];d=s['fermions'];alpha=sorted({l for c in contents(tuple(p['mu']),p['source_length']) for l in c});n=math.ceil((s['dimension']+16)/math.comb(d+2,d));return Frame(3,n,d+2,[d],alpha,seed=2026111917+7919*src,p=P)

def source_basis(p,words,src,a):
    from native_su_parallel_v1 import evaluate
    from su_modular_rowspace_v1 import RowSpace
    s=p['sources'][src];i=SOURCES.index(src);d=s['dimension'];fixed=[w for k,w in words if k==i];maps.setup(a.runtime.resolve(),p,list(TARGETS),list(SOURCES),P);frame=basis_frame(p,src)
    rng=random.Random(2026111923+src*1009);space=RowSpace(frame.ncols,P);sel=[];seen=set(fixed);fixed_input=list(fixed);attempts=0;started=time.monotonic()
    def add(block,require_all=False):
        expr=[maps.project(maps.algebra.parity_project({w:1},1),maps.charge,P) for w in block];v=evaluate(expr,frame,a.workers);ix=space.extend(v);assert not require_all or ix==list(range(len(block)));sel.extend(block[j] for j in ix)
    for group in batches(fixed_input,48):add(group,True)
    # Source point evaluations, independent of the differential evaluations, verify a complement.
    def checkpoint():
        fp=a.out/f'basis_{src}.json.gz';save_z(fp,dict(source=src,dimension=d,fixed=len(fixed),words=sel,columns=list(map(int,space.pivots)),frame_width=frame.ncols,seed=2026111917+7919*src,prime=P));save_json(a.out/f'basis_{src}.json',dict(source=src,dimension=d,rank=len(sel),fixed=len(fixed),complement=len(sel)-len(fixed),candidates_drawn=attempts,sha256=fsha(fp),complete=len(sel)==d,minor_sha256=INPUT));return len(sel)
    last=time.monotonic();checkpoint()
    while len(sel)<d and time.monotonic()-started<a.seconds:
        block=[]
        while len(block)<min(64, max(16, 2*(d-len(sel)))):
            w=maps.random_word(rng.choice(maps.choices[i]),rng,6);attempts+=1
            if w is not None and w not in seen:seen.add(w);block.append(w)
        add(block);assert len(sel)<=d
        if time.monotonic()-last>45:print(src,checkpoint(),d,flush=True);last=time.monotonic()
    print(src,checkpoint(),d,flush=True)

def init_pool(runtime,p,t):maps.setup(runtime,p,[t],list(SOURCES),P)

def side(p,words,t,a):
    active=[i for i,(src,_) in enumerate(words) if t in p['sources'][SOURCES[src]]['targets']];tasks=[words[i] for i in active];assert active
    init=(a.runtime.resolve(),p,t);init_pool(*init);width=maps.frames[0].ncols;base=a.out/f'side_{t}';base.mkdir(exist_ok=True)
    # A file consists only of rows with structural incidence, always in the original pivot order.
    manifest=dict(target=t,prime=P,plan_sha256=PLAN,minor_sha256=INPUT,pivot_rows=active,total_pivot_rows=len(words),width=width,batch=8)
    marker=base/'manifest.json'
    if marker.exists():assert json.loads(marker.read_bytes())==manifest
    else:save_json(marker,manifest)
    missing=[j for j in range(0,len(tasks),8) if not (base/(str(j)+'.npy')).exists()]
    with ProcessPoolExecutor(max_workers=a.workers,initializer=init_pool,initargs=init) as pool:
        # Bounded outstanding packets, even for the widest evaluation frames.
        for window in batches(missing,2*a.workers):
            vals=pool.map(maps.values,[tasks[j:j+8] for j in window]);
            for j,v in zip(window,vals):
                tmp=base/(str(j)+'.tmp');f=base/(str(j)+'.npy');assert v.shape==(min(8,len(tasks)-j),width)
                with tmp.open('wb') as wr:np.save(wr,pack(v),allow_pickle=False)
                tmp.replace(f)
    h=hashlib.sha256();nnz=0
    for j in range(0,len(tasks),8):
        row=np.load(base/(str(j)+'.npy'),allow_pickle=False);assert row.shape==(min(8,len(tasks)-j),width) and row.dtype==np.dtype('<u4') and int(row.max())<P
        h.update(row.tobytes());nnz+=np.count_nonzero(row)
    r=dict(target=t,prime=P,active_pivot_rows=len(active),width=width,entries=len(active)*width,nonzero=int(nnz),entries_le_u4_sha256=h.hexdigest(),manifest_sha256=fsha(marker),minor_sha256=INPUT);save_json(base/'result.json',r);print(json.dumps(r),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['factor','basis','side']);ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--plan',type=Path,default=maps.HERE/'classical_order1.json');ap.add_argument('--minor',type=Path,default=maps.HERE/'other_maps/task_5/minor.json.gz');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--profile',type=int);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--seconds',type=int,default=1800);a=ap.parse_args();assert a.workers>0;sys.path.insert(0,str(a.runtime.resolve()));a.out.mkdir(exist_ok=True,parents=True);p,words,cols=read_inputs(a.plan,a.minor)
    if a.mode=='factor':make_factor(p,words,cols,a)
    elif a.mode=='basis':assert a.profile in SOURCES;source_basis(p,words,a.profile,a)
    else:assert a.profile in OTHERS;side(p,words,a.profile,a)
