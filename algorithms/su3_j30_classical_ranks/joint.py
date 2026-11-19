"""Evaluate Q0 and Cartan restriction of the same ordered length-nine SU(3) seed.

The row indices are basis-recipe indices. Blocks use production target frames at
length ten, and traceless diagonal frames at length nine. SHA-256 pins precede
all unpickling. Only complete row intervals are committed to disk.
"""
import argparse
import hashlib
import json
import math
import pickle
import sys
from pathlib import Path
import numpy as np

P = 1_000_033
PLAN = '58a440c904dda13ed2f0a34ee3e26543edbd607049dd2774e8af6dadaef65ae6'
BASE = '3cd48eea41430f6cd08449f651f729db8125417ce04b1d6b8dbd6f4b77867149'
INPUTS = 'd261f2fc038af8fede3646c720b84d4d3e498cf6365b6c4990dd6e2b93fb26cc'
RT = 'runtime/work/sp3-native.nosync'
HERE = Path(__file__).resolve().parent

def sha(f):
    with Path(f).open('rb') as b: return hashlib.file_digest(b, 'sha256').hexdigest()

def digest(x):
    return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def save(f,x):
    tmp=f.with_suffix(f.suffix+'.part'); tmp.write_text(json.dumps(x,indent=2)+'\n'); tmp.replace(f)

def load_plan():
    f=HERE/'joint_plan.json'; assert __debug__ and sha(f)==PLAN
    return json.loads(f.read_bytes())

def inputs(root):
    # The original cluster staging manifest has extra job inputs. Compare the
    # actual complete runtime and basis subset with the published manifest pin.
    full=json.loads((root/'SOURCE_MANIFEST.json').read_bytes())
    names=sorted(n for n in full if (n.startswith('runtime/') and n.endswith(('.py','.cpp','.json'))) or n.startswith('second_basis/shard0/targets'))
    assert len(names)==221 and all(not Path(n).is_absolute() and '..' not in Path(n).parts for n in names)
    actual={n:sha(root/n) for n in names}; assert all(actual[n]==full[n] for n in names) and digest(actual)==INPUTS
    return INPUTS

def runtime(root):
    sys.path.insert(0,str(root.resolve()/RT)); import su3_second_correction_map_v1 as maps
    assert maps.PRIME==P; return maps

def recipes(root,s):
    base=root/'second_basis/shard0'/('targets%02d'%s['index']); f=base/'independent_seeds.pkl'
    assert sha(f)==s['recipe_sha256']; a=pickle.loads(f.read_bytes()); assert len(a)==s['dimension']; return a

def cartan_frame(maps,plan,i,points=0):
    s=plan['sources'][i]; d=s['fermions']; assert s['cartan_dimension']>0
    k=d+2; points=points or math.ceil((s['cartan_dimension']+32)/math.comb(k,d))
    _,letters=maps.alphabets(tuple(plan['mu']),9,10); seed=2026111919+7919*i
    f=maps.Frame(3,points,k,[d],letters,seed=seed,p=P,cartan=True)
    desc=dict(N=3,prime=P,k=k,degree=d,points=points,columns=f.ncols,seed=seed,letters=letters,cartan=True)
    return f,desc

def frames(maps,plan,i,points=0):
    s=plan['sources'][i]
    values={f'q{t:02d}':maps.target_frame(plan,t,P,points) for t in s['targets']}
    if s['cartan_dimension']:values['c']=cartan_frame(maps,plan,i,points)
    return values

def run_rows(maps,plan,i,words,fr,workers):
    from native_su_parallel_v1 import evaluate,evaluate_profile
    s=plan['sources'][i]; sector=maps.sector(tuple(plan['mu']),9)
    expr=[maps.project(maps.algebra.parity_project({w:1},-1),sector,P) for w in words]
    assert all(expr)
    def bounded(xs,f,evaluator):
        bound=max(1,8_000_000//f.ncols)
        if sum(map(len,xs))<=bound:return evaluator(xs)
        if len(xs)>1:
            mid=len(xs)//2;return np.vstack((bounded(xs[:mid],f,evaluator),bounded(xs[mid:],f,evaluator)))
        terms=list(xs[0].items());total=np.zeros((1,f.ncols),dtype=np.int64)
        for j in range(0,len(terms),bound): total=(total+evaluator([dict(terms[j:j+bound])]))%P
        return total
    out={}
    for key,(f,_) in fr.items():
        if key=='c': callback=lambda xs:evaluate(xs,f,workers)
        else:
            t=plan['targets'][int(key[1:])];callback=lambda xs:evaluate_profile(xs,f,s['profile'],t['profile'],workers)
        v=bounded(expr,f,callback);assert v.shape==(len(words),f.ncols) and v.min()>=0 and v.max()<P;out[key]=v.astype('<u4')
    return out

def make_guard(plan,i,fr):
    return dict(plan=PLAN,inputs=INPUTS,driver=sha(__file__),source=i,recipe=plan['sources'][i]['recipe_sha256'],frames={k:digest(v) for k,(_,v) in fr.items()})

def interval(root,out,plan,maps,i,lo,hi,workers,batch,cache):
    if i not in cache:
        fs=frames(maps,plan,i); cache.clear();cache[i]=(recipes(root,plan['sources'][i]),fs,make_guard(plan,i,fs))
    words,fs,guard=cache[i];assert 0<=lo<hi<=len(words) and 1<=batch<=32
    folder=out/f's{i:02d}'/f'{lo:05d}_{hi:05d}';folder.mkdir(parents=True,exist_ok=True);receipt=folder/'result.json'
    if receipt.exists():
        r=json.loads(receipt.read_bytes());assert r['guard']==guard and [r['start'],r['stop']]==[lo,hi] and set(r['matrices'])==set(fs)
        for k,x in r['matrices'].items(): assert sha(folder/(k+'.npy'))==x['sha256']
        return r
    arr={k:np.empty((hi-lo,f.ncols),dtype='<u4') for k,(f,_) in fs.items()}
    for a in range(lo,hi,batch):
        b=min(a+batch,hi);rows=run_rows(maps,plan,i,words[a:b],fs,workers)
        for k,v in rows.items():arr[k][a-lo:b-lo]=v
    r=dict(guard=guard,start=lo,stop=hi,matrices={})
    for k,v in arr.items():
        f=folder/(k+'.npy');tmp=f.with_suffix('.part')
        with tmp.open('wb') as buf:np.save(buf,v,allow_pickle=False)
        tmp.replace(f);r['matrices'][k]=dict(shape=list(v.shape),sha256=sha(f),nonzero=int(np.count_nonzero(v)))
    save(receipt,r);return r

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('root',type=Path);ap.add_argument('output',type=Path);ap.add_argument('source',type=int);ap.add_argument('start',type=int);ap.add_argument('stop',type=int);ap.add_argument('--workers',type=int,default=4);ap.add_argument('--batch',type=int,default=8);a=ap.parse_args();inputs(a.root);print(json.dumps(interval(a.root,a.output,load_plan(),runtime(a.root),a.source,a.start,a.stop,a.workers,a.batch,{})))
