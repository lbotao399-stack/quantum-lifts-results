"""Evaluate the two untouched components of -Q1 O0 on reproducible SU(3) frames."""
import argparse, gzip, hashlib, json, sys
from pathlib import Path
from collections import Counter
from functools import lru_cache
import numpy as np
import projected_q0 as maps
from coupled_first_replay import PRIME, RHS_SHA, coeff, make_frame, vector_sha
from coupled_kernel_fit import one_frame

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()))
    import sp_environment
    from sp_full_profiles import profile
    from g2_full_classical_screen_v1 import contents
    from native_su_values_v1 import evaluate
    root=a.root; a.out.mkdir(parents=True,exist_ok=True)
    plan=json.loads((root/'classical_order1.json').read_bytes())['incoming'];mu=tuple(plan['mu']);targets=(6,10)
    alphabet=sorted({l for k in (10,11) for c in contents(mu,k) for l in c}); cfg={t:plan['targets'][t] for t in targets};pats={maps.pattern(r):t for t,r in cfg.items()}
    assert all(set(plan['sources'][s]['targets']).isdisjoint(targets) for s in (19,20,21,22))
    frames={t:make_frame(r,alphabet,a.points,88304773+229*t) for t,r in cfg.items()};ds={t:np.zeros(f.ncols,dtype=np.int64) for t,f in frames.items()};buf={t:{} for t in targets};checked=Counter();ns=Counter();h=hashlib.sha256();total=0
    with gzip.open(root/'negative_Q1.jsonl.gz','rb') as src:
        for line in src:
            h.update(line);total+=1;c,v=json.loads(line);w=maps.decode(v);t=pats.get(profile(tuple(l for tr in w for l in tr)))
            if t is None:continue
            n=coeff(c);ns[t]+=1; assert w not in buf[t];buf[t][w]=n
            if checked[t]<3:
                f=one_frame(frames[t]);v=f.eval_words([w],dtype=np.int64)[0]*n%PRIME;assert np.array_equal(v,evaluate([{w:n}],f)[0]);checked[t]+=1
            if len(buf[t])>=256:
                ds[t]=(ds[t]+evaluate([buf[t]],frames[t])[0])%PRIME;buf[t]={}
    assert total==517468 and h.hexdigest()==RHS_SHA
    rec={}
    for t in targets:
        f=frames[t];arr=(ds[t]+evaluate([buf[t]],f)[0])%PRIME;assert ns[t]==plan['first_rhs_trace_terms'][str(t)] and checked[t]==3
        path=a.out/f'fresh_{t}.npy';np.save(path,arr.astype('<u4'),allow_pickle=False)
        nz=np.flatnonzero(arr);j=int(nz[0]) if len(nz) else None
        rec[str(t)]=dict(target=t,seed=88304773+229*t,prime=PRIME,noncommuting_points=f.P,exterior_coordinates_per_point=f.cpp,width=f.ncols,rhs_terms=ns[t],residual_nonzero=len(nz),first_nonzero=None if j is None else dict(index=j,value=int(arr[j])),vectors_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),residual_sha256=vector_sha(arr),rhs_stream_sha256=RHS_SHA,old_correction_has_no_path=True,python_exterior_controls=checked[t]);(a.out/f'fresh_{t}.json').write_text(json.dumps(rec[str(t)],indent=2)+'\n')
    print(json.dumps(rec),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=maps.HERE);p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--out',type=Path,required=True);p.add_argument('--points',type=int,default=12);run(p.parse_args())
