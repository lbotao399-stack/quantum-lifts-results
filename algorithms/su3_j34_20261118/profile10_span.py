"""Test the entire archived first-correction span on a level-34 profile-10 frame.

Rows are in expanded_edges.json.gz order; a zero frozen coefficient never filters a row.
Point blocks bound evaluation memory. A verified separating functional excludes the
selected field-valued directions in the parent problem. It does not exhaust the sources.
"""
import argparse, json, sys, os
from concurrent.futures import ProcessPoolExecutor
from copy import copy
from pathlib import Path
import numpy as np
import profile10_validation as v
from all_coordinate_rank import matrix, pivot_indices
from full_strand import fsha, harray, packed, putjson

P=v.PRIME

def part(frame,i,j):
    f=copy(frame);f.P=j-i;f.ncols=f.P*f.cpp
    f.L={l:(d,x[i:j]) for l,(d,x) in frame.L.items()}
    return f

def setup(rt,root,cert,points,seed):
    global plan,rows,coef,fr,charge,alg,project,derivative,evaluate,tpat,frames
    sys.path.insert(0,str(rt));import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from native_su_values_v1 import evaluate
    plan,rows,coef=v.inputs(root,cert);fr,charge=v.frame(plan,points,seed)
    tpat=v.maps.pattern(plan['targets'][10]);frames=[part(fr,i,min(i+16,points)) for i in range(0,points,16)]

def compute_batch(ix):
    s=rows[ix[0]][0];assert all(rows[i][0]==s for i in ix)
    pat=v.maps.pattern(plan['sources'][s]);es=[project(alg.parity_project({rows[i][1]:1},1),charge,P) for i in ix]
    vals=packed(np.concatenate([derivative(es,f,pat,tpat) for f in frames],axis=1))
    return ix,vals

def build(a):
    init=(a.runtime.resolve(),a.root,a.cert,a.points,a.seed);setup(*init)
    d=a.work; d.mkdir(parents=True,exist_ok=True); cfg=dict(seed=a.seed,points=a.points,prime=P,recipes=v.RECIPES_SHA,plan=v.PLAN_SHA,rank_record=v.RANK_SHA,rhs=v.RHS_SHA,program=fsha(Path(__file__)))
    manifest=d/'config.json'
    if manifest.exists():assert json.loads(manifest.read_bytes())==cfg
    else:putjson(manifest,cfg)
    rhs_path=d/'rhs.npy'
    if not rhs_path.exists():
        terms=v.rhs_terms(a.root,plan);b=np.zeros(fr.ncols,dtype=np.int64)
        for batch in v.batches(list(terms.items()),128):b=(b+np.concatenate([evaluate([dict(batch)],f)[0] for f in frames]))%P
        with (d/'rhs.tmp').open('wb') as out:np.save(out,packed(b),allow_pickle=False)
        (d/'rhs.tmp').replace(rhs_path);print('RHS complete',len(terms),flush=True)
    tasks=[]
    for s in v.EDGE:
        ix=[i for i,(k,w) in enumerate(rows) if k==s];tasks+=list(v.batches(ix,4))
    def filename(ix):return d/('rows_%04d.npy'%ix[0])
    todo=[q for q in tasks if not filename(q).exists()]
    with ProcessPoolExecutor(max_workers=a.workers,initializer=setup,initargs=init) as pool:
        for n,(ix,w) in enumerate(pool.map(compute_batch,todo)):
            assert w.shape==(len(ix),fr.ncols);dest=filename(ix);temp=dest.with_suffix('.tmp')
            with temp.open('wb') as out:np.save(out,w,allow_pickle=False)
            temp.replace(dest)
            if (n+1)%40==0:print('additional row batches',n+1,'of',len(todo),flush=True)
    b=np.load(rhs_path,allow_pickle=False);A=np.empty((len(rows),fr.ncols),dtype='<u4');used=np.zeros(len(rows),dtype=bool)
    for ix in tasks:
        u=np.load(filename(ix),allow_pickle=False);assert u.shape==(len(ix),fr.ncols) and u.dtype==A.dtype and not used[ix].any();used[ix]=True;A[ix]=u
    assert used.all() and b.shape==(fr.ncols,) and b.dtype==A.dtype and max(A.max(),b.max())<P
    checks={};from sp_full_profiles import profile
    for s in v.EDGE:
        i=next(i for i,(k,w) in enumerate(rows) if k==s and A[i,:fr.cpp].any());expr=project(alg.parity_project({rows[i][1]:1},1),charge,P);z={w:int(c)%P for w,c in alg.differential(expr).items() if c%P and profile(tuple(l for tr in w for l in tr))==tpat};f=part(fr,0,1)
        np.testing.assert_array_equal(A[i,:fr.cpp],evaluate([z],f)[0]);checks[str(s)]=i
    if (a.points,a.seed)==(3,761519013):
        old=json.loads((a.cert/'profile10_fresh.json').read_bytes());np.testing.assert_array_equal(b,old['rhs']);np.testing.assert_array_equal(coef@A.astype(np.int64)%P,old['lhs'])
    np.save(d/'matrix.npy',A,allow_pickle=False)
    info=dict(**cfg,shape=list(A.shape),rhs_data_sha256=harray(b),matrix_data_sha256=harray(A),direct_checks=checks,recipes_by_profile={str(s):sum(k==s for k,w in rows) for s in v.EDGE},other_structurally_zero_rows=1519)
    putjson(d/'frame.json',info);return info

def rank_test(d,out):
    info=json.loads((d/'frame.json').read_bytes());A=np.load(d/'matrix.npy',allow_pickle=False);b=np.load(d/'rhs.npy',allow_pickle=False)
    assert info['matrix_data_sha256']==harray(A) and info['rhs_data_sha256']==harray(b) and len(A)==2752 and info['shape']==list(A.shape)
    r=int(matrix(A).rank());u=int(matrix(np.vstack([A,b])).rank()); I=pivot_indices(matrix(A).transpose());sel=A[I];J=pivot_indices(matrix(sel));assert r==len(I)==len(J)
    square=matrix(sel[:,J]);det=int(square.det());assert det and r*(P-1)**2<2**63
    c=np.asarray([int(x[0]) for x in square.transpose().solve(matrix(b[J,None])).tolist()],dtype=np.int64);res=(b.astype(np.int64)-c@sel.astype(np.int64))%P;assert not res[J].any() and u==r+int(res.any())
    info.update(rank=r,augmented_rank=u,pivot_rows=I,pivot_columns=J,pivot_minor_determinant=det,membership=(r==u))
    if res.any():
        j=int(np.flatnonzero(res)[0]);dual=np.zeros(A.shape[1],dtype=np.int64);dual[j]=1;dual[J]=[-int(z[0]) for z in square.solve(matrix(sel[:,j,None])).tolist()];dual%=P
        assert all(not (block.astype(np.int64)@dual%P).any() for block in np.array_split(A,32));pair=int(b.astype(np.int64)@dual%P);assert pair==int(res[j])!=0
        info.update(annihilator=dual.tolist(),pairing=pair,extra_column=j,augmented_minor_determinant=det*pair%P,annihilated_recipe_rows=len(A))
    else:info.update(coefficients_on_pivot_rows=c.tolist())
    putjson(out,info);print('ranks',r,u,'pair',info.get('pairing'),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=v.maps.HERE);ap.add_argument('--cert',type=Path,default=v.maps.HERE/'allcoord');ap.add_argument('--runtime',type=Path,default=v.maps.RUNTIME);ap.add_argument('--work',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--points',type=int,default=280);ap.add_argument('--seed',type=int,default=419870023);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--rank-only',action='store_true');a=ap.parse_args();assert __debug__ and a.points>0 and a.workers>0
    if not a.rank_only:build(a)
    rank_test(a.work,a.output)
