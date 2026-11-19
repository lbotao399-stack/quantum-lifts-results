"""Re-evaluate SU(3) target recipes in the matching production frame modulo 1000033.

A selected minor uses original ordered seed rows followed, when required, by
identified packed Q0 image rows. Each extension retains the signed Schur minor.
In-memory computation is restarted from verified evaluation chunks after interruption.
"""
import argparse
import hashlib
import json
import pickle
import sys
from pathlib import Path
from zipfile import ZipFile

import numpy as np

P = 1000033
FIRST_ZIP = '866710743d70aa63cc3882b9ead5421e07dc0ae71d7d4585f0a3ec9aa511df72'
PLANS = {1: ('su3_correction_profiles_v1.json','b4eb0388bbaa1a6dc655e4b2d877ba48e771d00a89cc7a76d48f4278219570eb'),
         2: ('su3_second_correction_profiles_v1.json','b8c669fa244e851f38a25b91186c3091e98a939935242805b0f71a9f8c00d609')}

def digest(blob): return hashlib.sha256(blob).hexdigest()

def sha(path):
    with path.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def load(path): return json.loads(path.read_bytes())

def save(path,obj):
    temp=path.with_suffix('.json.part');temp.write_text(json.dumps(obj,separators=(',',':'))+'\n');temp.replace(path)

def field():
    from flint import nmod_mat
    return nmod_mat

def pivots(matrix,rank):
    pc=[];c=0
    for r in range(rank):
        while c<matrix.ncols() and not matrix[r,c]: c+=1
        assert c<matrix.ncols();pc.append(c);c+=1
    return pc

class Minor:
    def __init__(self,width):
        self.width=width;self.cols=[];self.rows=[];self.panels=[];self.det=1
        self.basis=np.empty((0,width),dtype=np.int64)
    def extend(self,rows,labels):
        from exact_blas_rows import matmul_mod
        fl=field(); rows=np.asarray(rows,dtype=np.int64)
        assert rows.shape==(len(labels),self.width) and np.all((rows>=0)&(rows<P))
        if not len(rows):return
        if self.cols: rows=(rows-matmul_mod(rows[:,self.cols],self.basis,P))%P
        if not rows.any():return
        rr,h=fl(rows.tolist(),P).rref();h=int(h);cc=pivots(rr,h)
        small=fl(rows[:,cc].T.tolist(),P);rr,h2=small.rref();assert int(h2)==h
        ii=pivots(rr,h);a=fl(rows[np.ix_(ii,cc)].tolist(),P);d=int(a.det());assert d
        new=matmul_mod(np.array(a.inv().tolist(),dtype=np.int64),rows[ii],P)
        assert np.array_equal(new[:,cc],np.eye(h,dtype=np.int64)) and not new[:,self.cols].any()
        if self.cols:self.basis=(self.basis-matmul_mod(self.basis[:,cc],new,P))%P
        self.basis=np.vstack((self.basis,new));self.cols.extend(cc);self.rows.extend(labels[i] for i in ii)
        self.panels.append([h,d]);self.det=self.det*d%P
    def result(self):
        return dict(rank=len(self.cols),rows=self.rows,columns=self.cols,panel_orders_and_determinants=self.panels,
                    determinant_in_selected_order=self.det)

def control():
    fl=field();rng=np.random.default_rng(72);checks=0
    for h,w in [(1,1),(8,6),(11,13),(16,12)]:
        for rk in range(0,min(h,w)+1):
            rows=(rng.integers(P,size=(h,rk))@rng.integers(P,size=(rk,w)))%P
            rows=rows[::-1,::-1].copy(); mat=fl(rows.tolist(),P)
            for b in (1,3,8):
                m=Minor(w)
                for i in range(0,h,b):m.extend(rows[i:i+b],list(range(i,min(i+b,h))))
                assert len(m.rows)==int(mat.rank())
                assert (int(fl(rows[np.ix_(m.rows,m.cols)].tolist(),P).det()) if m.cols else 1)==m.det
                checks+=1
    return checks

def evaluate_bounded(expressions,frame,workers):
    from native_su_parallel_v1 import evaluate
    # Bound the native word-value matrix without dropping a trace term.
    limit=max(1,8_000_000//frame.ncols)
    if sum(map(len,expressions))<=limit:return evaluate(expressions,frame,workers)
    if len(expressions)>1:
        k=len(expressions)//2
        return np.vstack((evaluate_bounded(expressions[:k],frame,workers),evaluate_bounded(expressions[k:],frame,workers)))
    terms=list(expressions[0].items());total=np.zeros((1,frame.ncols),dtype=np.int64)
    for lo in range(0,len(terms),limit):total=(total+evaluate([dict(terms[lo:lo+limit])],frame,workers))%P
    return total

def run(root,first_zip,out,equation,index,workers,batch):
    import su3_correction_map_v2 as first
    import su3_second_correction_map_v1 as second
    import su_trace_adapter_v1 as algebra
    from native_su_projection_v1 import project
    from g2_full_classical_screen_v1 import sector
    assert __debug__ and first.PRIME==second.PRIME==P and workers>=1 and 1<=batch<=64
    runtime=root/'runtime/work/sp3-native.nosync';nm,hash0=PLANS[equation];plan=load(runtime/nm);assert sha(runtime/nm)==hash0
    t=plan['targets'][index];assert t['index']==index and plan['parity']==-1
    length=10 if equation==1 else plan['target_length'];mu=tuple(plan['mu'])
    rhs=root/'output'/('first' if equation==1 else 'second_rhs')/'rhs'/('profile%02d'%index)
    meta=load(rhs/'manifest.json');stamp=load(rhs/'result.json');assert meta['target']==t
    frame,desc=first.target_frame(t) if equation==1 else second.target_frame(plan,index)
    desc=json.loads(json.dumps(desc));assert desc==meta['frame'] and desc['prime']==P and desc['columns']==frame.ncols==stamp['columns']
    assert stamp['status']==f'complete_negative_Q{equation}_evaluations' and sha(rhs/'rhs.npy')==stamp['sha256']
    if equation==1:
        assert sha(first_zip)==FIRST_ZIP
        with ZipFile(first_zip) as z: raw=z.read(f'first{index:02}.pkl')
    else:
        folder=root/'second_basis/shard0'/('targets%02d'%index);record=load(folder/'result.json');base=load(folder/'manifest.json')
        assert base['task']==t and base['mu']==list(mu) and base['length']==length and base['role']=='targets' and base['parity']==-1
        raw=(folder/'independent_seeds.pkl').read_bytes();assert digest(raw)==record['seeds_sha256'] and record['dimension']==t['dimension']
    seeds=pickle.loads(raw);assert len(seeds)<=t['dimension'] and (equation==1 or len(seeds)==t['dimension'])
    dst=out/('e%d_t%02d'%(equation,index));dst.mkdir(parents=True,exist_ok=True);q=dst/'manifest.json'
    manifest=dict(equation=equation,target=t,prime=P,frame=desc,recipe_sha256=digest(raw),recipe_count=len(seeds),batch=batch,
                  plan_sha256=hash0,rhs_manifest_sha256=sha(rhs/'manifest.json'),rhs_sha256=stamp['sha256'],driver_sha256=sha(Path(__file__)),
                  evaluator_sha256=sha(runtime/'native_su_parallel_v1.py'),projector_sha256=sha(runtime/'native_su_projection_v1.py'))
    if q.exists():assert load(q)==manifest
    else:save(q,manifest)
    m=Minor(frame.ncols);chunks=[];maps=[]
    for lo in range(0,len(seeds),batch):
        hi=min(len(seeds),lo+batch);name=f'rows{lo:05}_{hi:05}';fp=dst/(name+'.npy');ck=dst/(name+'.json')
        if not ck.exists():
            expr=[project(algebra.parity_project({word:1},-1),sector(mu,length),P) for word in seeds[lo:hi]]
            vals=evaluate_bounded(expr,frame,workers);assert vals.shape==(hi-lo,frame.ncols) and np.all((vals>=0)&(vals<P))
            tmp=fp.with_suffix('.part.npy');np.save(tmp,vals.astype('<u4'));tmp.replace(fp);save(ck,dict(start=lo,stop=hi,sha256=sha(fp)))
        rec=load(ck);assert (rec['start'],rec['stop'],rec['sha256'])==(lo,hi,sha(fp));vals=np.load(fp,allow_pickle=False);assert vals.dtype==np.dtype('<u4')
        m.extend(vals,[['seed',i] for i in range(lo,hi)]);chunks.append(rec)
        save(dst/'progress.json',dict(equation=equation,target=index,evaluated=hi,seed_rank=len(m.cols),dimension=t['dimension']))
    seed_rank=len(m.cols)
    if equation==1 and seed_rank<t['dimension']:
        for s in plan['sources']:
            if index not in s['targets'] or len(m.cols)==t['dimension']:continue
            n=s['dimension'];i=s['index'];folder=root/'output/compacted/first/maps'/('profile%02d'%i)/('rows%05d_%05d'%(0,n));meta=load(folder/'manifest.json');assert meta['source']==s and meta['targets'][str(index)]==desc and meta['prime']==P
            prefix='t%02d_%05d_%05d'%(index,0,n);fp=folder/(prefix+'.npy');rec=load(folder/(prefix+'.json'));assert (rec['source'],rec['target'],rec['start'],rec['stop'],rec['columns'],rec['sha256'])==(i,index,0,n,frame.ncols,sha(fp))
            v=np.load(fp,mmap_mode='r',allow_pickle=False);assert v.shape==(n,frame.ncols) and v.dtype==np.dtype('<u4')
            for lo in range(0,n,batch):
                hi=min(lo+batch,n);m.extend(v[lo:hi],[['image',i,j] for j in range(lo,hi)])
                if len(m.cols)==t['dimension']:break
            maps.append(dict(source=i,sha256=rec['sha256'],rows_used_through=hi))
    assert len(m.cols)<=t['dimension'] and m.det
    result=dict(equation=equation,target=index,prime=P,dimension=t['dimension'],width=frame.ncols,injective=len(m.cols)==t['dimension'],seed_rank=seed_rank,
                witness=m.result(),manifest_sha256=sha(q),chunks=chunks,packed_Q0_sources=maps)
    save(dst/'result.json',result);print(json.dumps({k:result[k] for k in ('equation','target','dimension','injective','seed_rank')}|{'rank':len(m.cols),'determinant':m.det}),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--first-zip',type=Path)
    ap.add_argument('--out',type=Path);ap.add_argument('--equation',type=int,choices=(1,2));ap.add_argument('--index',type=int);ap.add_argument('--workers',type=int,default=8);ap.add_argument('--batch',type=int,default=16);ap.add_argument('--control',action='store_true');args=ap.parse_args()
    sys.path.insert(0,str(args.root/'runtime/work/sp3-native.nosync'));import sp_environment
    if args.control:print(json.dumps({'signed_minor_controls':control()}))
    else:run(args.root,args.first_zip,args.out,args.equation,args.index,args.workers,args.batch)
