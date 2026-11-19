"""Minors at p=1000033 for the recorded Cartan frames and length-eight predecessors.

The Cartan rows are taken from the original Q0/Cartan run. Predecessor rows
are freshly evaluated from the pinned recipes at p (old basis receipts used
another prime). No random rank estimate or numerical tolerance is used.
"""
import argparse
import hashlib
import json
import math
import pickle
import sys
from pathlib import Path
import numpy as np

P = 1000033
PLAN = '58a440c904dda13ed2f0a34ee3e26543edbd607049dd2774e8af6dadaef65ae6'
SOURCE = 'b8c669fa244e851f38a25b91186c3091e98a939935242805b0f71a9f8c00d609'
JREC = '8d9f8f13fc8b5e640d2cb48a269d4ad65a997667adbd6267d4c15b9b3789e677'
LISTS = 'e7c39b5545afd3c57a2104a6e451d24a2908127aa3c3136681620983564a26a4'

def sha(p):
    with Path(p).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def read(p): return json.loads(p.read_bytes())
def signed(w): return math.prod(x[1] for x in w['panel_orders_and_determinants']) % P
def save(p,data):
    t=p.with_suffix('.part');t.write_text(json.dumps(data,separators=(',',':'))+'\n');t.replace(p)

def run(science,method,stamps,base,out,kind,i,workers=4):
    rt=science/'runtime/work/sp3-native.nosync';sys.path.insert(0,str(rt));import sp_environment
    import su3_second_correction_map_v1 as maps
    from reproject import Minor,evaluate_bounded
    from flint import ctx
    ctx.threads=workers
    from joint import digest,load_plan,cartan_frame
    assert __debug__ and 1<=workers<=16 and kind in ('c','p')
    plan=load_plan(); assert sha(method/'joint_plan.json')==PLAN and sha(rt/'su3_second_correction_profiles_v1.json')==SOURCE
    old=read(rt/'su3_second_correction_profiles_v1.json');mu=tuple(plan['mu']);assert mu==tuple(old['mu'])
    key=f'{kind}{i:02d}';out.mkdir(parents=True,exist_ok=True);prov={};base=Path(base)
    if kind=='c':
        assert sha(stamps/'joint_receipt.json')==JREC
        src=plan['sources'][i];dim=src['cartan_dimension'];frame,desc=cartan_frame(maps,plan,i);m=Minor(frame.ncols)
        jf=base/'su3_j30_joint_rank_20261119'/'output'/(key+'.json');stage=next(s for s in read(stamps/'joint_receipt.json')['stages'] if s['stage']==key)
        assert sha(jf)==stage['stage_sha256'];g=read(jf)['guard'];assert g['key']==key and g['sources']==[i] and g['width']==frame.ncols and g['plan']==PLAN
        recs=[];lo0=0; used=[]
        for f in sorted((base/'su3_j30_joint_map_20261119'/'output'/f's{i:02d}').glob('*_*/result.json')):
            r=read(f);b=r['matrices']['c'];lo,hi=r['start'],r['stop'];assert lo0==lo<hi<=src['dimension'] and r['guard']['recipe']==src['recipe_sha256'] and r['guard']['frames']['c']==digest(desc) and r['guard']['plan']==PLAN
            lo0=hi;recs.append([i,lo,hi,b['sha256']]);assert b['shape']==[hi-lo,frame.ncols]
            if len(m.cols)<dim:
                path=f.parent/'c.npy';assert sha(path)==b['sha256'];a=np.load(path,allow_pickle=False);assert a.shape==tuple(b['shape']) and a.dtype==np.dtype('<u4');m.extend(a,list(range(lo,hi)));used.append([lo,hi,b['sha256']])
        assert lo0==src['dimension'] and digest(recs)==g['blocks'];prov=dict(recipe_sha256=src['recipe_sha256'],stage_sha256=sha(jf),used_chunks=used)
    else:
        assert sha(stamps/'source_order.json')==LISTS
        s=old['sources'][i];row=next(r for r in read(stamps/'source_order.json')['profiles'] if (r['equation'],r['index'])==(2,i));dim=s['dimension'];length=old['source_length'];d=s['fermions'];f=science/'second_basis/shard0'/f'sources{i:02d}'/'independent_seeds.pkl'
        assert row['dimension']==dim and row['profile']==s['profile'] and sha(f)==row['seeds_sha256'] and length==8
        words=pickle.loads(f.read_bytes());assert len(words)==dim
        alphabet=sorted({L for c in maps.contents(mu,length) for L in c});k=d+2;points=math.ceil((dim+32)/math.comb(k,d));seed=2026111723+i*7919
        frame=maps.Frame(3,points,k,[d],alphabet,seed=seed,p=P);desc=dict(N=3,k=k,points=points,degree=d,alphabet=alphabet,seed=seed,prime=P);m=Minor(frame.ncols);sector=maps.sector(mu,length)
        for lo in range(0,dim,32):
            hi=min(dim,lo+32);expr=[maps.project(maps.algebra.parity_project({w:1},-1),sector,P) for w in words[lo:hi]]
            a=evaluate_bounded(expr,frame,workers);m.extend(a,list(range(lo,hi)))
        prov=dict(recipe_sha256=sha(f))
    w=m.result();assert 0<len(m.cols)<=dim and signed(w)==m.det!=0 and len(m.rows)==len(set(m.rows))==len(set(m.cols))==len(m.cols)
    rec=dict(key=key,prime=P,dimension=dim,width=frame.ncols,rank=w['rank'],saturated=w['rank']==dim,frame_sha256=digest(desc),witness=w,**prov)
    save(out/(key+'.json'),rec);print(key,rec['rank'],dim,flush=True);return rec

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('science',type=Path);a.add_argument('method',type=Path);a.add_argument('stamps',type=Path);a.add_argument('base',type=Path);a.add_argument('out',type=Path);a.add_argument('kind',choices=['c','p']);a.add_argument('index',type=int);a.add_argument('--workers',type=int,default=4);x=a.parse_args();sys.path.insert(0,str(x.method.resolve()));sys.path.insert(0,str((x.method.parent/'su3_j30_frame_audit').resolve()));run(x.science.resolve(),x.method.resolve(),x.stamps,x.base,x.out,x.kind,x.index,x.workers)
