"""Fit the four incoming families of target 10, using one coefficient vector at all 17 targets.

Keep the last four of the twelve independently generated points of each target out of
elimination. The six recorded nonzero right-hand sides are residuals after the old
profile-16 correction. Incidence zeros are exact; point tests never prove a cochain identity.
"""
import argparse, gzip, hashlib, json, random, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from flint import nmod_mat
import projected_q0 as maps
from coupled_first_replay import PRIME as P, PLAN_SHA, RHS_SHA, make_frame, vector_sha
from coupled_kernel_fit import one_frame

SOURCE=(7,9,16,17)
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def put(path,v): path.write_text(json.dumps(v,indent=2)+'\n')
def arr_sha(a): return vector_sha(a)

def run(root,runtime,out):
    sys.path.insert(0,str(runtime));import su_trace_adapter_v1 as algebra
    from g2_full_classical_screen_v1 import contents,sector
    from sp_full_profiles import profile
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as block
    from native_su_values_v1 import evaluate
    from su_modular_rowspace_v1 import RowSpace
    from su3_stream_rank_v3 import random_word
    q=root/'classical_order1.json';assert sha(q)==PLAN_SHA
    p=json.loads(q.read_bytes())['incoming'];mu=tuple(p['mu']);charge=sector(mu,10)
    assert tuple(r['index'] for r in p['sources'] if 10 in r['targets'])==SOURCE
    ts=sorted({t for s in SOURCE for t in p['sources'][s]['targets']});assert ts==[2,3,6,8,10,11,13,14]
    alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});cs=defaultdict(list)
    for c in contents(mu,10): cs[profile(c)].append(c)
    out.mkdir(parents=True,exist_ok=True);frames={};rhs=[];train=[];hold=[];offset=0;bounds={};old={};checks=[];old_hashes={}
    for t,r in enumerate(p['targets']):
        f=make_frame(r,alpha,12,88304773+229*t);frames[t]=f
        if t in (6,10):
            d=root/'joint/edges';v=np.load(d/f'fresh_{t}.npy',allow_pickle=False);rec=json.loads((d/f'fresh_{t}.json').read_bytes());assert rec['rhs_stream_sha256']==RHS_SHA
        elif t in (13,14,15,16):
            d=root/'first_profile16/coupled';v=np.load(d/f'fresh_{t}.npy',allow_pickle=False)[2];rec=json.loads((d/f'fresh_{t}.json').read_bytes());assert rec['formal_rhs_stream_sha256']==RHS_SHA
        else:
            assert str(t) not in p['first_rhs_trace_terms'] and all(t not in p['sources'][s]['targets'] for s in (19,20,21,22));v=np.zeros(f.ncols,dtype=np.int64);rec=None
        if rec:
            assert rec['seed']==88304773+229*t and rec['prime']==P and rec['width']==f.ncols and rec['residual_sha256']==arr_sha(v) and rec['vectors_sha256']==sha(d/f'fresh_{t}.npy');old_hashes[str(t)]=rec['vectors_sha256']
        old[t]=np.asarray(v,dtype=np.int64);rhs.extend(v);bounds[t]=(offset,offset+f.ncols)
        if t in ts:
            train.extend(range(offset,offset+8*f.cpp));hold.extend(range(offset+8*f.cpp,offset+f.ncols))
        offset+=f.ncols
    rhs=np.array(rhs,dtype=np.int64);d=len(train);space=RowSpace(d,P);M=[];selected=[];seen=set();batch=0;controls=set();rng=random.Random(3407101617)
    while len(selected)<d and batch<700:
        s=SOURCE[batch%len(SOURCE)];pat=maps.pattern(p['sources'][s]);words=[]
        while len(words)<12:
            w=random_word(rng.choice(cs[pat]),rng,6)
            if w is not None and (s,w) not in seen:words.append(w);seen.add((s,w))
        expr=[project(algebra.parity_project({w:1},1),charge,P) for w in words];V=np.zeros((len(words),len(rhs)),dtype=np.int64)
        for t in p['sources'][s]['targets']:
            f=frames[t];lo,hi=bounds[t];tar=maps.pattern(p['targets'][t]);V[:,lo:hi]=block(expr,f,pat,tar)
            if (s,t) not in controls:
                for i,e in enumerate(expr):
                    if not np.any(V[i,lo:lo+f.cpp]):continue
                    u={w:c%P for w,c in algebra.differential(e).items() if c%P and profile(tuple(l for tr in w for l in tr))==tar};np.testing.assert_array_equal(V[i,lo:lo+f.cpp],evaluate([u],one_frame(f))[0]);checks.append([s,t,arr_sha(V[i,lo:lo+f.cpp])]);controls.add((s,t));break
        keep=space.extend(V[:,train]);selected.extend((s,words[i]) for i in keep);M.extend(V[i].astype('<u4') for i in keep);batch+=1
        if batch%32==0: print('batches',batch,'rank',len(selected),'of',d,flush=True)
        if len(controls)==sum(len(p['sources'][s]['targets']) for s in SOURCE) and len(selected)==d:break
    assert len(controls)==sum(len(p['sources'][s]['targets']) for s in SOURCE)
    M=np.asarray(M,dtype='<u4');r=len(selected);cols=[train[i] for i in space.pivots];minor=nmod_mat(M[:,cols].tolist(),P);det=int(minor.det());assert det
    c=np.array([int(z[0]) for z in minor.transpose().solve(nmod_mat([[int(rhs[i])] for i in cols],P)).tolist()],dtype=np.int64);assert r*(P-1)**2<2**63
    change=c@M.astype(np.int64)%P;res=(rhs-change)%P;assert not np.any(res[cols]);stats={}
    for t,f in frames.items():
        lo,hi=bounds[t];n=8*f.cpp;v=res[lo:hi];stats[str(t)]=dict(reachable=t in ts,fit_width=n,withheld_width=4*f.cpp,before=int(np.count_nonzero(rhs[lo:hi])),training_nonzero=int(np.count_nonzero(v[:n])),withheld_nonzero=int(np.count_nonzero(v[n:])),residual_sha256=arr_sha(v));
        if t not in ts:np.testing.assert_array_equal(v,old[t]);assert not np.any(change[lo:hi])
    recipe=out/'recipes.json.gz';recipe.write_bytes(gzip.compress(json.dumps(selected,separators=(',',':')).encode(),mtime=0));arr=out/'matrix.npy';np.save(arr,M,allow_pickle=False);sol=out/'solution.npz';np.savez_compressed(sol,coefficients=c.astype('<u4'),rhs=rhs.astype('<u4'),residual=res.astype('<u4'),columns=np.array(cols,dtype='<u4'))
    result=dict(prime=P,plan_sha256=PLAN_SHA,baseline_hashes=old_hashes,sources=SOURCE,reachable=ts,rank=r,training_width=d,withheld_width=len(hold),whole_width=len(rhs),minor_determinant=det,minor_sha256=arr_sha(M[:,cols]),training_nonzero=int(np.count_nonzero(res[train])),withheld_nonzero=int(np.count_nonzero(res[hold])),all_target_stats=stats,input_recipes_by_source=dict(Counter(s for s,w in selected)),nonzero_coefficients_by_source=dict(Counter(s for (s,w),v in zip(selected,c) if v)),recipe_sha256=sha(recipe),matrix_sha256=sha(arr),solution_sha256=sha(sol),explicit_word_controls=checks,batches=batch,fit_points_per_target=8,withheld_points_per_target=4);put(out/'result.json',result);print({k:result[k] for k in ('rank','training_width','minor_determinant','training_nonzero','withheld_nonzero')},flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=maps.HERE);ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();run(a.root,a.runtime.resolve(),a.out)
