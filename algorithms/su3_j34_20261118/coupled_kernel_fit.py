"""Fit the saved residual on a strict subset of independent coordinates using coupled partial boundaries.

Every column changes the classical residual only on profiles 13,14,15. Replaying
on withheld coordinates distinguishes a partial fit from a joint solution.
"""
import argparse, gzip, hashlib, json, random, sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from copy import copy
from functools import lru_cache
from pathlib import Path
import numpy as np
from flint import nmod_mat
from kernel_strands import HERE, FAMILIES, describe, load
from projected_q0 import decode, pattern, RUNTIME

P=1000003
TARGETS=(13,14,15)

def sha(data):return hashlib.sha256(data).hexdigest()

@lru_cache(maxsize=16384)
def word_profile(w):return profile(tuple(x for trace in w for x in trace))

def init(root,runtime):
    global algebra,project,contents,sector,profile,evaluate,eval_block,random_word,frames,choices,geom,pre,plan,sids,tids,charge
    sys.path.insert(0,str(runtime));import su_trace_adapter_v1 as algebra
    from native_su_projection_v1 import project
    from g2_full_classical_screen_v1 import contents,sector
    from sp_full_profiles import profile
    from native_su_values_v1 import evaluate
    from native_su_profile_derivation_v1 import evaluate as eval_block
    from su3_stream_rank_v3 import random_word
    from coupled_first_replay import make_frame
    word_profile.cache_clear(); geom=describe(root)['families'];pre,plan=load(root);mu=tuple(pre['mu']);charge=sector(mu,9)
    sids={pattern(r):r['index'] for r in plan['sources']};tids={pattern(r):r['index'] for r in plan['targets']};c9=defaultdict(list)
    for c in contents(mu,9):c9[profile(c)].append(c)
    choices={i:c9[pattern(pre['sources'][i])] for i in FAMILIES}
    alphabet=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});frames={}
    for t in TARGETS:
        rec=json.loads((root/'first_profile16/coupled'/f'fresh_{t}.json').read_bytes());assert (rec['mode'],rec['target'],rec['width'],rec['prime'])==('fresh',t,252,P)
        frames[t]=make_frame(plan['targets'][t],alphabet,rec['noncommuting_points'],rec['seed'])

def keep_parts(blocks, expression, ids):
    res={i:{} for i in blocks}
    for w,c in expression.items():
        i=ids.get(word_profile(w));
        if i in res and c%P:res[i][w]=int(c)%P
    return res

def trace_split(w,i):
    rec=geom[str(i)]; image=keep_parts(rec['kept'],algebra.differential({w:1}),sids)
    return {s:project(e,sector(tuple(plan['mu']),10),P) if e else {} for s,e in image.items()}

def one_frame(f):
    first=copy(f);first.P=1;first.ncols=f.cpp;first.L={l:(d,v[:1]) for l,(d,v) in f.L.items()}
    for k in tuple(first.__dict__):
        if k.startswith('_native_su_'):delattr(first,k)
    return first

def projected_parts_after_two(i,w):
    # Independent order: highest-weight project only after BOTH unprojected differentials.
    rec=geom[str(i)];parity=algebra.parity_project({w:1},1);pieces=keep_parts(rec['kept'],algebra.differential(parity),sids)
    complete=Counter();
    for e in pieces.values():complete.update(algebra.differential(e))
    grouped=keep_parts(TARGETS,complete,tids)
    assert not any(c%P and tids.get(word_profile(v)) not in rec['shared_final_targets'] for v,c in complete.items())
    return {t:project(e,sector(tuple(plan['mu']),11),P) if e else {} for t,e in grouped.items()}

def rows(packet,control=False):
    if not packet:return np.empty((0,756),dtype=np.int64)
    i=packet[0][0];assert all(j==i for j,_ in packet);rec=geom[str(i)]
    mids=[]
    for _,w in packet:
        # Apply charge conjugation once. Both gradings, both projectors, and Q0 are equivariant.
        par=algebra.parity_project({w:1},1); mid=keep_parts(rec['kept'],algebra.differential(par),sids)
        mids.append({s:project(e,sector(tuple(plan['mu']),10),P) if e else {} for s,e in mid.items()})
    out=np.zeros((len(packet),len(TARGETS),252),dtype=np.int64)
    for j,t in enumerate(TARGETS):
        if t not in rec['shared_final_targets']:continue
        for s in rec['kept']:
            if t in plan['sources'][s]['targets']:
                out[:,j]=(out[:,j]+eval_block([d[s] for d in mids],frames[t],pattern(plan['sources'][s]),pattern(plan['targets'][t])))%P
    if control:
        # Formal nilpotence before projection, integer signed traces (no finite-rank evaluations).
        w=packet[0][1];assert not algebra.differential(algebra.differential({w:1}))
        test=projected_parts_after_two(i,w)
        for j,t in enumerate(TARGETS):
            f=one_frame(frames[t]);assert np.array_equal(evaluate([test[t]],f)[0],out[0,j,:f.cpp]),(i,t)
    return out.reshape(len(packet),-1)

def sample(batch):
    i=FAMILIES[batch%len(FAMILIES)];rng=random.Random(908731+104729*batch); sample=[]
    while len(sample)<6:
        w=random_word(rng.choice(choices[i]),rng,6)
        if w is not None:sample.append((i,w))
    return sample,rows(sample,control=batch<len(FAMILIES))

def saved_rows(recipes):
    # Preserve original row order across profile groups.
    data=np.zeros((len(recipes),756),dtype=np.int64)
    for i in FAMILIES:
        ix=[k for k,(idx,_) in enumerate(recipes) if idx==i]
        for j in range(0,len(ix),8):
            keep=ix[j:j+8];data[keep]=rows([recipes[k] for k in keep],control=(j==0))
    return data

def sha_file(p):return sha(p.read_bytes())

def run(a):
    from contextlib import nullcontext
    root=a.root.resolve();runtime=a.runtime.resolve();init(root,runtime)
    from su_modular_rowspace_v1 import RowSpace
    inputs={};rhs=[]
    for t in TARGETS:
        f=root/'first_profile16/coupled'/f'fresh_{t}.npy';j=json.loads(f.with_suffix('.json').read_bytes());inputs[str(t)]=sha_file(f);assert inputs[str(t)]==j['vectors_sha256'];v=np.load(f,allow_pickle=False).astype(np.int64);assert np.array_equal(v[2],(v[0]-v[1])%P);rhs.extend(v[2]);
    rhs=np.array(rhs,dtype=np.int64);assert a.fit_points in range(1,13) and a.workers>=1 and a.batches>=7
    cols=[k for k in range(756) if (k%252)//21<a.fit_points];other=[k for k in range(756) if k not in cols];out=a.out;out.mkdir(parents=True,exist_ok=True)
    if a.replay:
        witness=json.loads(gzip.decompress(a.replay.read_bytes()));assert inputs==witness['residual_array_sha256'] and witness['columns_to_fit']==cols and witness['prime']==P
        recipes=[(i,decode(w)) for i,w in witness['recipes']];basis=saved_rows(recipes);piv=witness['pivot_columns'];weights=np.array(witness['coefficients'],dtype=np.int64);n=len(piv)
    else:
        space=RowSpace(len(cols),P);recipes=[];memo=[];end=0
        pool=ProcessPoolExecutor(max_workers=a.workers,initializer=init,initargs=(root,runtime)) if a.workers>1 else None
        with pool if pool else nullcontext():
            while end<a.batches and space.rows.shape[0]<len(cols):
                span=range(end,min(end+max(7,2*a.workers),a.batches));packets=list(pool.map(sample,span) if pool else map(sample,span));keys=[w for k,_ in packets for w in k];vectors=np.vstack([v for _,v in packets]);ind=space.extend(vectors[:,cols]);recipes.extend(keys[k] for k in ind);memo.extend(vectors[ind]);end=span.stop
                print('batches',end,'rank',len(recipes),'of',len(cols),flush=True)
        if not recipes:raise ValueError('no nonzero correction')
        basis=np.stack(memo);piv=[cols[i] for i in space.pivots];n=len(piv)
        weights=np.array(nmod_mat(basis[:,piv].T.tolist(),P).solve(nmod_mat(n,1,rhs[piv].tolist(),P)).tolist(),dtype=np.int64).ravel()
    assert len(set(recipes))==len(weights)==n and list(sorted(set(piv)))==sorted(piv)
    determinant=int(nmod_mat(basis[:,piv].tolist(),P).det());assert determinant
    assert len(weights)*(P-1)**2<2**63;after=(rhs-weights@basis)%P;assert not after[piv].any()
    counts={str(t):dict(fit_nonzero=int(np.count_nonzero(after[[k for k in cols if k//252==j]])),held_out_nonzero=int(np.count_nonzero(after[[k for k in other if k//252==j]]))) for j,t in enumerate(TARGETS)}
    facts=dict(prime=P,rank=n,determinant=determinant,fit_coordinates=len(cols),withheld_coordinates=len(other),before_nonzero=int(np.count_nonzero(rhs)),fit_nonzero=int(np.count_nonzero(after[cols])),withheld_nonzero=int(np.count_nonzero(after[other])),by_target=counts,recipes_by_predecessor=dict(sorted(Counter(str(s) for s,_ in recipes).items())),nonzero_weights=int(np.count_nonzero(weights)),result_vector_sha256=sha(after.astype('<u4').tobytes()),residual_array_sha256=inputs)
    if a.replay:
        for k,v in facts.items():assert v==witness[k],k
        print(json.dumps(facts),flush=True);return
    cert=dict(**facts,columns_to_fit=cols,pivot_columns=piv,recipes=recipes,coefficients=weights.tolist());file=out/'witness.json.gz';file.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));(out/'result.json').write_text(json.dumps(dict(**facts,witness_sha256=sha_file(file)),indent=2)+'\n');print(json.dumps(facts),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=HERE);p.add_argument('--runtime',type=Path,default=RUNTIME);p.add_argument('--workers',type=int,default=1);p.add_argument('--fit-points',type=int,default=2);p.add_argument('--batches',type=int,default=140);p.add_argument('--out',type=Path,required=True);p.add_argument('--replay',type=Path);run(p.parse_args())
