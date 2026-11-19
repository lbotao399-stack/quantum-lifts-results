"""Complete target 11 from range-checked source-10 chunks. Never read a partial array as input."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import schur_side_images as side
sc, m, maps = side.sc, side.m, side.maps
T, S, N, W, STEP, PARTS = 11, 10, 1554, 35500, 195, 8


def inputs(a):
    sys.path.insert(0,str(a.runtime.resolve()))
    p, words, _=m.read_inputs(a.plan,a.minor); maps.setup(a.runtime.resolve(),p,[T],list(m.SOURCES),m.P)
    f=maps.frames[0]; assert f.ncols==W
    active,pv,pvh=side.pivot_images(a,T,words,p,W)
    assert len(active)==413 and pvh=='48533deeaa1d145e8c47ca53fa3614a696e309d28f204b5503947aa4def0ba34'
    records=[]; B=None
    for s in m.SOURCES:
        idx,rows,k,_,basis=sc.source_basis(a,p,words,s); free=rows[k:]
        cp=a.coords/f'coefficients_{s}.npy'; cert=json.loads((a.coords/f'coordinates_{s}.json').read_bytes());cs=m.fsha(cp)
        assert cert['coefficient_file_sha256']==cs and cert['source_basis_file_sha256']==basis and cert['source_minor_determinant'] and cert['prime']==m.P and cert['nonzero_residuals']==0 and cert['pivot_minor_sha256']==m.INPUT
        c=np.load(cp,allow_pickle=False);assert c.shape==(len(free),len(words)) and c.dtype==np.dtype('<u4');co=np.ascontiguousarray(c[:,active]);use=T in p['sources'][s]['targets'];nn=int(np.count_nonzero(co));zero=not len(free) or (not nn and not use)
        assert zero==(s!=S);r=dict(source=s,free=len(free),width=W,source_basis_sha256=basis,coefficients_sha256=cs,free_reaches_target=use,active_pivot_coefficients=nn,structural_zero=zero)
        if zero:r['nonzero_entries']=0
        else:
            assert (len(free),nn)==(N,528719);B=(idx,free,co)
        records.append(r)
    return p,words,f,active,pv,pvh,records,B


def direct(p,f,idx,words,co,pv,row,point):
    sub=side.cut(f,point,point+1); sl=slice(point*f.cpp,(point+1)*f.cpp)
    v=sc.grouped_action(p,sub,T,(idx,words[row]),maps.charge)
    corr=np.array([sum(int(a)*int(b) for a,b in zip(co[row],col))%m.P for col in pv[:,sl].T],dtype=np.int64)
    return (v-corr)%m.P


def chunk(a):
    from exact_blas_rows import matmul_mod
    from native_su_profile_derivation_v1 import evaluate
    p,_,f,_,pv,pvh,records,(idx,words,co)=inputs(a);start=a.part*STEP;end=min(start+STEP,N);assert 0<=start<end<=N
    out=a.work;out.mkdir(exist_ok=True,parents=True);dest=out/f'part{a.part}.npy';rp=dest.with_suffix('.json')
    if rp.is_file():
        r=json.loads(rp.read_bytes());assert (r['interval'],r['sha256'],r['coefficients_sha256'],r['pivot_hash'],r['plan'])==([start,end],m.fsha(dest),records[-1]['coefficients_sha256'],pvh,m.PLAN);return
    tmp=dest.with_suffix('.tmp.npy');x=np.lib.format.open_memmap(tmp,mode='w+',dtype='<u4',shape=(end-start,W));points=[(i,min(i+256,f.P)) for i in range(0,f.P,256)];frames=[side.cut(f,i,j) for i,j in points];source,target=maps.pattern(p['sources'][S]),maps.pattern(p['targets'][T]);count=0
    for i in range(start,end,12):
        j=min(i+12,end);expr=[maps.project(maps.algebra.parity_project({w:1},1),maps.charge,m.P) for w in words[i:j]];v=-matmul_mod(co[i:j],pv,m.P)
        for fr,(lo,hi) in zip(frames,points):v[:,lo*f.cpp:hi*f.cpp]+=evaluate(expr,fr,source,target)
        v%=m.P;x[i-start:j-start]=v;count+=int(np.count_nonzero(v))
    x.flush();del x;disk=np.load(tmp,mmap_mode='r',allow_pickle=False);checks=[]
    for row,point in [(start,0),((start+end)//2,f.P//2),(end-1,f.P-1)]:
        np.testing.assert_array_equal(disk[row-start,point*f.cpp:(point+1)*f.cpp],direct(p,f,idx,words,co,pv,row,point));checks.append([row,point,f.cpp])
    assert int(np.count_nonzero(disk))==count;del disk;tmp.replace(dest);r=dict(interval=[start,end],prime=m.P,plan=m.PLAN,source_basis_sha256=records[-1]['source_basis_sha256'],coefficients_sha256=records[-1]['coefficients_sha256'],pivot_hash=pvh,sha256=m.fsha(dest),nonzero=count,direct=checks);m.save_json(rp,r);print(json.dumps(r),flush=True)


def merge(a):
    p,allwords,f,active,pv,pvh,records,(idx,words,co)=inputs(a);out=a.out/f'target_{T}';out.mkdir(exist_ok=True,parents=True);dest=out/f'source_{S}.npy';tmp=a.work/'assembled.tmp.npy';x=np.lib.format.open_memmap(tmp,mode='w+',dtype='<u4',shape=(N,W));end=0;count=0;manifest=[]
    for part in range(PARTS):
        fp=a.work/f'part{part}.npy';r=json.loads(fp.with_suffix('.json').read_bytes());s,j=r['interval'];assert s==end==part*STEP and j==min(N,s+STEP) and (r['plan'],r['prime'],r['coefficients_sha256'],r['source_basis_sha256'],r['pivot_hash'])==(m.PLAN,m.P,records[-1]['coefficients_sha256'],records[-1]['source_basis_sha256'],pvh) and r['sha256']==m.fsha(fp);y=np.load(fp,mmap_mode='r',allow_pickle=False);assert y.shape==(j-s,W) and y.dtype==np.dtype('<u4') and np.count_nonzero(y)==r['nonzero'];x[s:j]=y;count+=r['nonzero'];end=j;manifest.append(r)
    assert end==N;x.flush();del x;y=np.load(tmp,mmap_mode='r',allow_pickle=False);assert np.count_nonzero(y)==count;point=f.P//2;checks=[0,N//2,N-1];sl=slice(point*f.cpp,(point+1)*f.cpp)
    for j in checks:np.testing.assert_array_equal(y[j,sl],direct(p,f,idx,words,co,pv,j,point))
    del y;sha=m.fsha(tmp);assert not dest.exists() or m.fsha(dest)==sha;tmp.replace(dest);records[-1].update(path=dest.name,sha256=sha,nonzero_entries=count,direct_rows=checks,direct_point=point,direct_coordinates=len(checks)*f.cpp);tests=[0,len(active)//2,len(active)-1];one=side.cut(f,point,point+1)
    for j in tests:np.testing.assert_array_equal(pv[j,sl],sc.grouped_action(p,one,T,allwords[active[j]],maps.charge))
    result=dict(target=T,prime=m.P,minor_sha256=m.INPUT,plan_sha256=m.PLAN,pivot_entries_sha256=pvh,pivot_direct_rows=[active[j] for j in tests],point=point,width=W,rows=records,implicit_or_stored_entries=sum(r['free']*W for r in records),sign='free image minus full pivot combination');m.save_json(a.work/'chunks.json',manifest);final=out/'result.json';assert not final.exists() or json.loads(final.read_bytes())==result;m.save_json(final,result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['chunk','merge']);ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--part',type=int,choices=range(PARTS));ap.add_argument('--plan',type=Path,default=maps.HERE/'classical_order1.json');ap.add_argument('--minor',type=Path,default=maps.HERE/'other_maps/task_5/minor.json.gz')
    for key in ('bases','coords','sides','work','out'):ap.add_argument('--'+key,type=Path,required=True)
    a=ap.parse_args();sys.path.insert(0,str(a.runtime.resolve()));assert __debug__;chunk(a) if a.mode=='chunk' else merge(a)
