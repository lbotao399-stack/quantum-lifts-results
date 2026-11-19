"""Prove injectivity of the length-ten evaluation frames in the R-007 sector.

The signed minor uses original Q0 images and, only if needed, independently
projected target words. Removing pivot columns between Schur updates keeps
large matrix products rectangular. Labels always refer to original rows.
"""
import argparse, json, math, random, sys
from pathlib import Path
import numpy as np
import joint

P=joint.P
JREC='8d9f8f13fc8b5e640d2cb48a269d4ad65a997667adbd6267d4c15b9b3789e677'

def read(p):return json.loads(p.read_bytes())

def take_pivots(a,h):
    j=0;out=[]
    for i in range(h):
        while not a[i,j]:j+=1
        out.append(j);j+=1
    return out

class Panels:
    def __init__(self,w):
        self.live=list(range(w));self.ops=[];self.ii=[];self.jj=[];self.dets=[];self.det=1
    @staticmethod
    def update(a,cols,keep,norm):
        # 128 summands, all products and every partial sum are exact binary64 integers.
        if not len(a):return np.zeros((0,len(keep)),dtype=np.uint32)
        assert 0<len(cols)<=128 and len(cols)*(P-1)**2<2**53
        prod=np.ascontiguousarray(a[:,cols],dtype=np.float64)@norm
        return ((a[:,keep].astype(np.int64)-prod.astype(np.int64))%P).astype(np.uint32)
    def panel(self,a,labels):
        from flint import nmod_mat
        assert 0<len(a)<=128 and a.shape==(len(labels),len(self.live))
        if not a.any():return
        mat=nmod_mat(a.tolist(),P);red,h=mat.rref();h=int(h);jj=take_pivots(red,h)
        rr,h2=nmod_mat(a[:,jj].T.tolist(),P).rref();assert int(h2)==h;ii=take_pivots(rr,h)
        sq=nmod_mat(a[np.ix_(ii,jj)].tolist(),P);d=int(sq.det());assert d
        skip=set(jj);keep=[j for j in range(len(self.live)) if j not in skip];u=np.asarray([[int(x) for x in row] for row in sq.inv().tolist()],dtype=np.float64)
        norm=(u@np.ascontiguousarray(a[np.ix_(ii,keep)],dtype=np.float64)).astype(np.int64)%P
        op=(np.array(jj,dtype=np.intp),np.array(keep,dtype=np.intp),np.ascontiguousarray(norm,dtype=np.float64));assert not self.update(a,*op).any()
        self.ii.extend(labels[i] for i in ii);self.jj.extend(self.live[j] for j in jj);self.dets.append([h,d]);self.det=self.det*d%P;self.ops.append(op);self.live=[self.live[j] for j in keep];return op
    def fresh(self,a,labels):
        for op in self.ops:
            a=self.update(a,*op)
            if not a.any():return
        return self.panel(a,labels)
    def stored(self,A,labels,cap,panel=96,tile=512):
        for i in range(0,len(A),panel):
            if len(self.ii)==cap:break
            n=len(self.live);stop=min(i+panel,len(A));op=self.panel(A[i:stop,:n],labels[i:stop]);assert len(self.ii)<=cap
            if op is not None and len(self.ii)<cap:
                for lo in range(stop,len(A),tile):
                    z=A[lo:lo+tile,:n];z[:,:len(self.live)]=self.update(z,*op)
    def witness(self):return dict(rank=len(self.ii),rows=self.ii,columns=self.jj,panel_orders_and_determinants=self.dets,determinant_in_selected_order=self.det)

def control():
    from flint import nmod_mat
    checks=0;g=np.random.default_rng(251119)
    for shape,h in [((9,15),1),((24,27),14),((19,13),10),((17,19),9),((8,7),0)]:
        n,w=shape;A=g.integers(P,size=(n,h),dtype=np.int64)@g.integers(P,size=(h,w),dtype=np.int64)%P;A[0]=0;A[:,0]=0;A[2]=A[1];A=A.astype(np.uint32)
        for k in (1,4,7):
            m=Panels(w);a=A[:n//2].copy();m.stored(a,list(range(len(a))),w,panel=k,tile=3)
            for i in range(len(a),n,k):m.fresh(A[i:i+k],list(range(i,min(i+k,n))))
            assert len(m.ii)==int(nmod_mat(A.tolist(),P).rank());assert m.det==(int(nmod_mat(A[np.ix_(m.ii,m.jj)].tolist(),P).det()) if m.ii else 1);checks+=1
    return checks

def run(base,out,i,workers,attempt_cap):
    from reproject import evaluate_bounded
    root=base/'su3_prime1000033_v1_20261118';rt=root/joint.RT;sys.path.insert(0,str(rt));import sp_environment
    from flint import ctx
    ctx.threads=workers;from su3_stream_rank_v3 import random_word
    from sp_full_profiles import profile
    assert __debug__ and control()==15 and 1<=workers<=64
    plan=joint.load_plan();maps=joint.runtime(root);key=f'q{i:02d}';target=plan['targets'][i];dim=target['dimension'];assert target['index']==i
    frame,desc=maps.target_frame(plan,i,P);w=frame.ncols;idx=[s['index'] for s in plan['sources'] if i in s['targets']]
    old=base/'su3_j30_joint_rank_20261119/output';q=old/(key+'.json');assert joint.sha(old/'result.json')==JREC
    stage=next(x for x in read(old/'result.json')['stages'] if x['stage']==key);assert joint.sha(q)==stage['stage_sha256'];g=read(q)['guard'];assert g['key']==key and g['sources']==idx and g['width']==w and g['plan']==joint.PLAN and g['worker']==joint.sha(joint.__file__)
    inc=stage['rank_increment'];assert inc<dim and inc<=w
    out.mkdir(parents=True,exist_ok=True);raw=base/'su3_j30_joint_map_20261119/output';M=sum(plan['sources'][s]['dimension'] for s in idx);A=np.empty((M,w),dtype='<u4');labels=[];links=[];pos=0
    for si in idx:
        src=plan['sources'][si];at=0
        for f in sorted((raw/f's{si:02d}').glob('*_*/result.json')):
            meta=read(f);lo,hi=meta['start'],meta['stop'];v=meta['matrices'][key];fp=f.parent/(key+'.npy');guard=meta['guard'];assert at==lo<hi<=src['dimension'] and guard['recipe']==src['recipe_sha256'] and guard['frames'][key]==joint.digest(desc) and guard['plan']==joint.PLAN
            assert v['shape']==[hi-lo,w] and joint.sha(fp)==v['sha256'];rows=np.load(fp,allow_pickle=False,mmap_mode='r');assert rows.shape==(hi-lo,w) and rows.dtype==np.dtype('<u4') and rows.max()<P
            A[pos:pos+hi-lo]=rows;pos+=hi-lo;at=hi;labels.extend([si,j] for j in range(lo,hi));links.append([si,lo,hi,v['sha256']]);del rows
        assert at==src['dimension']
    assert pos==M==len(labels) and joint.digest(links)==g['blocks'];m=Panels(w);m.stored(A,labels,dim);del A
    image_rank=len(m.ii);assert inc<=image_rank<=dim
    rng=random.Random(2026112041+i*1009);mu=tuple(plan['mu']);sel=tuple((tuple(x),n) for x,n in target['profile']);choices=[c for c in maps.contents(mu,10) if profile(c)==sel];assert choices
    seen=set();pending=[];words=[];attempts=0;seed=2026112041+i*1009
    def store():
        joint.save(out/(key+'.json'),dict(key=key,prime=P,plan=joint.PLAN,stage_sha256=stage['stage_sha256'],blocks=g['blocks'],sources=idx,image_rows=M,frame_sha256=joint.digest(desc),width=w,dimension=dim,increment=inc,image_rank=image_rank,rank=len(m.ii),saturated=len(m.ii)==dim,attempts=attempts,sampler_seed=seed,witness=m.witness()))
    store();print(key,'image',image_rank,'of',dim,flush=True)
    for attempts in range(1,attempt_cap+1) if image_rank<dim else ():
        word=random_word(rng.choice(choices),rng,6)
        if word is not None and word not in seen:
            seen.add(word);expr=maps.project(maps.algebra.parity_project({word:1},-1),maps.sector(mu,10),P)
            if expr:words.append(['trace',word]);pending.append(expr)
        if len(pending)>=32 or (pending and attempts==attempt_cap):
            a=evaluate_bounded(pending,frame,workers);assert a.shape==(len(pending),w) and a.min()>=0 and a.max()<P;m.fresh(a.astype(np.uint32),words);assert len(m.ii)<=dim;pending=[];words=[]
            if len(m.ii)==dim:break
        if attempts%2048==0:store();print(key,'attempt',attempts,'rank',len(m.ii),flush=True)
    store();print(key,'complete',len(m.ii),dim,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('base',type=Path);ap.add_argument('out',type=Path);ap.add_argument('index',type=int);ap.add_argument('--workers',type=int,default=16);ap.add_argument('--attempts',type=int,default=300000);a=ap.parse_args();run(a.base,a.out,a.index,a.workers,a.attempts)
