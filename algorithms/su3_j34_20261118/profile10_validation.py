"""Test the frozen 12-point fit on profile 10. Also exports a bounded projected membership test.

The archived baseline and all 1,519 kernel rows have no source profile mapping to 10.
Only the four raw edge profiles can contribute. All words and coefficients precede the fresh seed.
"""
import argparse, gzip, hashlib, json, sys
from collections import Counter
from pathlib import Path
import numpy as np
import projected_q0 as maps
from coupled_first_replay import coeff, make_frame, PRIME, PLAN_SHA, RHS_SHA
from full_strand import fsha, harray, packed, putjson

RANK_SHA='2cfb5767c9fbf3600dd8f5f0c51fdc57de5b6e14a55e6d836d8e56cd9d1ec9e9'
RECIPES_SHA='04b7d746722b478aee58c3ed14cace7a04379842db44831c15f09a1df055aef0'
EDGE=(7,9,16,17)

def inputs(root,cert):
    plan=root/'classical_order1.json';assert fsha(plan)==PLAN_SHA
    p=json.loads(plan.read_bytes())['incoming'];assert [s['index'] for s in p['sources'] if 10 in s['targets']]==list(EDGE)
    # First Schur seeds and their pivot directions; extra seed; carried strands; fixed baseline.
    kernel={1,2,4,10,6,13,14,15,18};base={19,20,21,22};assert set(EDGE).isdisjoint(kernel|base)
    for s in kernel|base:assert 10 not in p['sources'][s]['targets']
    path=cert/'expanded_edges.json.gz'; assert fsha(path)==RECIPES_SHA
    recipes=[(s,maps.decode(w)) for s,w in json.loads(gzip.decompress(path.read_bytes()))];assert len(recipes)==len(set(recipes))==2752 and set(s for s,w in recipes)==set(EDGE)
    path=cert/'all_rank.json';assert fsha(path)==RANK_SHA;r=json.loads(path.read_bytes());assert r['shape']==[4271,2244] and r['rank']==r['augmented_rank']==1308 and r['membership']
    c=np.zeros(4271,dtype=np.int64);c[r['pivot_rows']]=r['coefficients_on_pivot_rows'];assert np.count_nonzero(c)==1308 and np.count_nonzero(c[:2752])==1056
    return p,recipes,c[:2752]

def rhs_terms(root,p):
    from sp_full_profiles import profile
    pat=maps.pattern(p['targets'][10]);h=hashlib.sha256();n=0;terms={}
    with gzip.open(root/'negative_Q1.jsonl.gz','rb') as f:
        for line in f:
            h.update(line);n+=1;c,w=json.loads(line);w=maps.decode(w)
            if profile(tuple(l for tr in w for l in tr))==pat:
                assert w not in terms;terms[w]=coeff(c)
    assert n==517468 and h.hexdigest()==RHS_SHA and len(terms)==p['first_rhs_trace_terms']['10']==27780
    return terms

def frame(p,points,seed):
    from g2_full_classical_screen_v1 import contents,sector
    mu=tuple(p['mu']);alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});return make_frame(p['targets'][10],alpha,points,seed),sector(mu,10)

def batches(a,n):
    for i in range(0,len(a),n):yield a[i:i+n]

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from native_su_values_v1 import evaluate
    from sp_full_profiles import profile
    p,ww,c=inputs(a.root,a.cert);f,charge=frame(p,a.points,a.seed);tpat=maps.pattern(p['targets'][10]);terms=rhs_terms(a.root,p);arr=np.zeros(f.ncols,dtype=np.int64)
    for b in batches(list(terms.items()),256):arr=(arr+evaluate([dict(b)],f)[0])%PRIME
    lhs=np.zeros(f.ncols,dtype=np.int64);stats={};checks={}
    for s in EDGE:
        ix=[i for i,(k,w) in enumerate(ww) if k==s and c[i]];stats[s]=len(ix);pat=maps.pattern(p['sources'][s]);hits=0
        for batch in batches(ix,16):
            ex=[project(alg.parity_project({ww[i][1]:1},1),charge,PRIME) for i in batch];v=derivative(ex,f,pat,tpat);lhs=(lhs+c[batch]@v)%PRIME
            if not hits:
                for k,e in enumerate(ex):
                    if not v[k].any():continue
                    z={w:int(j)%PRIME for w,j in alg.differential(e).items() if j%PRIME and profile(tuple(l for tr in w for l in tr))==tpat};np.testing.assert_array_equal(v[k],evaluate([z],f)[0]);hits+=1;checks[s]=batch[k];break
        assert hits==1
    residual=(arr-lhs)%PRIME;bad=np.flatnonzero(residual);original=a.seed==88304773+229*10 and a.points==12
    if original:
        old=np.load(a.root/'joint/edges/fresh_10.npy',allow_pickle=False);np.testing.assert_array_equal(old,arr);assert not len(bad)
    result=dict(prime=PRIME,seed=a.seed,points=a.points,profile=10,original_frame=original,rank_certificate_sha256=RANK_SHA,edge_recipes_sha256=RECIPES_SHA,formal_rhs_sha256=RHS_SHA,trace_terms=len(terms),active_edges=stats,structurally_absent_kernel_rows=1519,explicit_differential_checks=checks,coordinate_count=f.ncols,rhs=arr.tolist(),lhs=lhs.tolist(),residual=residual.tolist(),residual_nonzero=len(bad),first_index=int(bad[0]) if len(bad) else None,residual_sha256=harray(residual));putjson(a.out,result);print(a.seed,f.ncols,len(bad),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=maps.HERE);ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--cert',type=Path,default=maps.HERE/'allcoord');ap.add_argument('--seed',type=int,required=True);ap.add_argument('--points',type=int,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert __debug__ and a.points>0;run(a)
