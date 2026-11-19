"""Carry source-18 pivot coefficients from target 12 into target 11 on the same source rows."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import full_strand as fs
import schur_side_images as side
from shared_first import SIZES


def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));assert fs.fsha(a.plan)==fs.PLAN;p=json.loads(a.plan.read_bytes())['incoming'];src,t=18,11;q=12
    assert set(p['sources'][src]['targets'])=={11,12,14} and all(11 not in p['sources'][s]['targets'] for s in (13,14,15))
    b=a.full/'target_12';f=b/'result.json';r=json.loads(f.read_bytes());assert r['prime']==fs.P and r['strand_source']==src and r['target']==q and r['projected_rows']==3900
    wsrc,words=fs.readwords(a.base/'minor.json.gz',fs.STRANDS[q][2]);assert r['minor_sha256']==fs.STRANDS[q][2] and wsrc['map']['sources']==[src] and r['full_images']['recipes_sha256']==fs.harray_words(words)
    cp=b/'coefficients.npy';assert fs.fsha(cp)==r['coefficient_file_sha256'];C=np.load(cp,mmap_mode='r',allow_pickle=False);assert C.shape==(3900,len(words))
    d=a.reduced/'target_11';prev=json.loads((d/'result.json').read_bytes());assert (prev['prime'],prev['plan_sha256'],prev['target'],prev['width'])==(fs.P,fs.PLAN,t,35500)
    out=a.out/'target_11';out.mkdir(parents=True,exist_ok=True);x,info=fs.build(out/'source18_images.npy',words,p,t,src,a.runtime,a.workers);assert x.shape==(3799,35500)
    frame=fs.maps.frames[0];pts=[0,frame.P//2,frame.P-1];one=[side.cut(frame,k,k+1) for k in pts];rows=[];offset=0
    for old,(s,n) in zip(prev['rows'],SIZES[:4]):
        assert (old['source'],old['free'],old['width'])==(s,n,x.shape[1]);v=C[offset:offset+n];offset+=n;zero=not np.any(v);raw=None
        if not old['structural_zero']:
            fp=d/old['path'];assert fs.fsha(fp)==old['sha256'];raw=np.load(fp,mmap_mode='r',allow_pickle=False);assert raw.shape==(n,x.shape[1])
        structural=n==0 or (raw is None and zero);z=dict(source=s,free=n,width=x.shape[1],structural_zero=structural,old_structural_zero=old['structural_zero'],active_strand_coefficients=int(np.count_nonzero(v)))
        if structural: z['nonzero_entries']=0;rows.append(z);continue
        dst=out/f'source_{s}.npy';tmp=dst.with_suffix('.tmp.npy');a0=np.lib.format.open_memmap(tmp,mode='w+',dtype='<u4',shape=(n,x.shape[1]));total=0
        for i in range(0,n,64):
            h=min(64,n-i);oldv=0 if raw is None else raw[i:i+h].astype(np.int64);v0=(oldv-(0 if zero else fs.product(v[i:i+h],x)))%fs.P;a0[i:i+h]=v0;total+=int(np.count_nonzero(v0))
        a0.flush();del a0;y=np.load(tmp,mmap_mode='r',allow_pickle=False);checks=[]
        for row,pt in zip((0,n//2,n-1),pts):
            sl=slice(pt*frame.cpp,(pt+1)*frame.cpp);res=np.asarray([sum(int(c)*int(e) for c,e in zip(v[row],col))%fs.P for col in x[:,sl].T],dtype=np.int64);orig=0 if raw is None else raw[row,sl].astype(np.int64);np.testing.assert_array_equal(y[row,sl],(orig-res)%fs.P);checks.append([row,pt,frame.cpp])
        assert int(np.count_nonzero(y))==total;del y;tmp.replace(dst);z.update(path=dst.name,sha256=fs.fsha(dst),nonzero_entries=total,direct=checks);rows.append(z)
    assert offset==3900;outrec=dict(target=t,prime=fs.P,plan_sha256=fs.PLAN,width=x.shape[1],source18_carry_included=True,source18_target12_receipt_sha256=fs.fsha(f),initial_target11_receipt_sha256=fs.fsha(d/'result.json'),source18_coefficients_sha256=fs.fsha(cp),source18_images=info,rows=rows,sign='first Schur side minus the target-12 strand coefficients times the source-18 side on target 11');fs.putjson(out/'result.json',outrec);print(json.dumps(outrec),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--runtime',type=Path,default=fs.maps.RUNTIME);ap.add_argument('--plan',type=Path,default=fs.maps.HERE/'classical_order1.json');ap.add_argument('--workers',type=int,default=4)
    for k in ('base','full','out','reduced'):ap.add_argument('--'+k,type=Path,required=True)
    a=ap.parse_args();assert 1<=a.workers<=8 and __debug__;run(a)
