"""Incoming Q0 of the pinned SU(3) descendants. Each row is one projected predecessor."""
import argparse, hashlib, json, pickle, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
P = 1000033
PLANS = {'Tp21':'86f03d87b9534642cb01c6b65f081cc45b4d33a1c5b0a34308ffca96f8b013d4',
         'S20':'ec5856736aa6a2510ba4dd94473bf5952055b63c1832cec44395e60fc980a9bd'}
SELECTED = {'Tp21':(3,5,20), 'S20':(3,6,10,12,13)}
RECEIPT = '3366b625abc25263450155b5e9693ac2efe618989e36ccd3079e4788626edf52'
PINS = {'su3_profile_basis_general_v1.py':'28bbb36cd20d56679f2249af05d565d66376b8bd6b8054ed0cb87424ca9017e9',
        'su3_second_correction_map_v1.py':'4ae02802084ce5fe640c6f833b7440ffe2a240c4938af0c76967cae4a893c1de',
        'native_su_profile_derivation_v1.py':'ad7b72a703e6ed181b00e4a0b54bba76462f18cf459d5ae3b5d044c8cefd175a',
        'native_su_projection_v1.py':'8be0e7be2f1fa86948e2d06c40fe830f895370a215d245562f62ce449c0f8406'}

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f, 'sha256').hexdigest()
def get(f):return json.loads(f.read_bytes())
def save(f,x):
    f.parent.mkdir(parents=True,exist_ok=True);t=f.with_suffix(f.suffix+'.part'); t.write_text(json.dumps(x,indent=2)+'\n');t.replace(f)
def array(f,a):
    assert np.issubdtype(a.dtype,np.integer) and np.all((0<=a)&(a<P));t=f.with_suffix('.part')
    with t.open('wb') as out:np.save(out,a.astype('<u4'),allow_pickle=False)
    t.replace(f);return sha(f)
def runtime(rt):
    assert __debug__ and all(sha(rt/k)==v for k,v in PINS.items());sys.path.insert(0,str(rt.resolve()))
    import su3_profile_basis_general_v1 as b, su3_second_correction_map_v1 as m
    assert sha(b.__file__)==PINS['su3_profile_basis_general_v1.py'] and m.PRIME==P
    return b,m
def plan(op):
    f=HERE/(op+'_plan.json');assert sha(f)==PLANS[op];a=get(f)
    assert [a[k] for k in ('prime','source_length','target_length','parity','op')]==[P,8,9,-1,op]
    for j in range(len(a['targets'])):assert a['targets'][j]['sources']==[s['index'] for s in a['sources'] if j in s['targets']]
    return a,f
def frame(m,a,j):
    f,d=m.target_frame(a,j,a['prime']); h=hashlib.sha256(json.dumps(d,sort_keys=True,separators=(',',':')).encode())
    for l,(parity,v) in sorted(f.L.items()):
        h.update(bytes((l[0],l[1],sum(1<<i for i in l[2]),parity)));h.update(np.ascontiguousarray(v,dtype='<i8').tobytes())
    assert f.p==d['prime']==P;return f,h.hexdigest()
def rhs(rt,raw,out,receipt,op,j,workers):
    from inputs import records,unpack
    from receipt import filename
    from native_su_parallel_v1 import evaluate
    from fortuity.full.evalg import nferm
    _,m=rt;a,pl=plan(op);assert sha(receipt)==RECEIPT;r=get(receipt);d=r['descendants'][op];assert r['prime']==P and d['plan_sha256']==PLANS[op]
    target=a['targets'][j];s=d['profile_blocks'][j];file=raw/op/filename(target);assert s['target']==target['index']==j and sha(file)==s['sha256'] and file.stat().st_size==s['bytes'];f,h=frame(m,a,j)
    v=np.zeros((1,f.ncols),dtype=np.int64);part={};terms=0
    for k,c in records(file):
        w=unpack(k);assert nferm(w)==target['fermions'] and w not in part;part[w]=c;terms+=1
        if len(part)==512:v=(v+evaluate([part],f,workers))%P;part.clear()
    if part:v=(v+evaluate([part],f,workers))%P
    assert terms==s['terms'];out.mkdir(parents=True,exist_ok=True);fileout=out/f'{op}_t{j:02d}.npy';digest=array(fileout,v);save(fileout.with_suffix('.json'),dict(op=op,prime=P,target=j,frame=h,plan_sha256=PLANS[op],input_sha256=s['sha256'],terms=terms,columns=f.ncols,nonzero=int(np.count_nonzero(v)),values_sha256=digest));return terms

def basis(b,a,pl,role,i,out,workers,seconds):
    # The frozen driver's default is 1,000,003. Pin its process-local binding to the checked plan.
    b.PRIME=a['prime']; assert b.PRIME==P;folder=out/f'{role}{i:02d}';task=a[role][i]
    if not (folder/'result.json').is_file() or get(folder/'result.json')['rank_lower_bound']!=task['dimension']:
        b.run(i,folder,pl,role,seconds=seconds,batch=12,workers=workers)
    mf=get(folder/'manifest.json');r=get(folder/'result.json');assert mf['prime']==P and mf['plan_sha256']==sha(pl) and mf['task']==task and mf['parity']==-1 and mf['role']==role and mf['driver_sha256']==PINS['su3_profile_basis_general_v1.py']
    assert r['status']=='complete_rational_profile_basis_recipes' and r['rank_lower_bound']==r['dimension']==task['dimension'],(role,i,r['rank_lower_bound'],task['dimension']);f=folder/'independent_seeds.pkl';assert r['seeds_sha256']==sha(f)
    seeds=pickle.loads(f.read_bytes());assert len(seeds)==task['dimension'];return seeds,r['seeds_sha256']

def worker(rt,op,role,i,out,workers,seconds,selected=None):
    b,m=rt;a,pl=plan(op);out=out/op;seeds,hseed=basis(b,a,pl,role,i,out,workers,seconds);folder=out/f'{role}{i:02d}';task=a[role][i]
    if role=='targets':
        f,h=frame(m,a,i);mf=get(folder/'manifest.json');assert mf['columns']==f.ncols and mf['frame']['seed']==2026111723+i*7919
        # Compare every target input matrix before and after adding derivative letters.
        g=b.Frame(3,mf['frame']['points'],mf['frame']['k'],mf['frame']['degrees'],sorted({l for c in b.contents(tuple(a['mu']),9) for l in c}),seed=mf['frame']['seed'],p=P)
        assert all(f.L[l][0]==p and np.array_equal(f.L[l][1],v) for l,(p,v) in g.L.items());save(folder/'done.json',dict(op=op,role=role,prime=P,index=i,plan_sha256=PLANS[op],dimension=len(seeds),seeds_sha256=hseed,frame=h));return
    ids=SELECTED[op] if selected is None else tuple(selected)
    assert ids and len(set(ids))==len(ids) and all(0<=j<len(a['targets']) for j in ids)
    targets=sorted(set(task['targets'])&set(ids));assert targets
    frames={j:frame(m,a,j) for j in targets};blocks={j:np.empty((len(seeds),f.ncols),dtype='<u4') for j,(f,_) in frames.items()}; sector=b.sector(tuple(a['mu']),8);step=6
    from native_su_parallel_v1 import evaluate_profile
    for lo in range(0,len(seeds),step):
        xx=[b.project(b.algebra.parity_project({w:1},-1),sector,P) for w in seeds[lo:lo+step]];assert all(xx)
        for j,(f,_) in frames.items():
            y=evaluate_profile(xx,f,task['profile'],a['targets'][j]['profile'],workers);assert np.all((0<=y)&(y<P));blocks[j][lo:lo+len(xx)]=y
    meta={j:dict(frame=frames[j][1],columns=v.shape[1],sha256=array(folder/f't{j:02d}.npy',v)) for j,v in blocks.items()}
    save(folder/'done.json',dict(op=op,role=role,prime=P,index=i,plan_sha256=PLANS[op],dimension=len(seeds),seeds_sha256=hseed,targets=meta));print(op,role,i,len(seeds),targets,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('runtime',type=Path);p.add_argument('op',choices=PLANS);p.add_argument('role',choices=['rhs','sources','targets']);p.add_argument('index',type=int);p.add_argument('out',type=Path);p.add_argument('--raw',type=Path);p.add_argument('--receipt',type=Path);p.add_argument('--workers',type=int,default=8);p.add_argument('--seconds',type=int,default=1500);p.add_argument('--targets',type=int,nargs='+',help='Coupled target profile indices; defaults reproduce the original selection.');x=p.parse_args();rt=runtime(x.runtime)
    if x.role=='rhs':print(x.op,x.index,rhs(rt,x.raw,x.out,x.receipt,x.op,x.index,x.workers),flush=True)
    else:worker(rt,x.op,x.role,x.index,x.out,x.workers,x.seconds,x.targets)
