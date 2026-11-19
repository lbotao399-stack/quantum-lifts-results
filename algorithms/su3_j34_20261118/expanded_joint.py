"""Replay every side of the 1,519 selected kernel directions and extend the four target-10 sources.
The saved 12-point frames are distinct from the full frames that selected the kernel.
Only the first eight points enter selection or fitting; the final four are validation.
"""
import argparse, gzip, json, math, random, sys
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from flint import nmod_mat
import projected_q0 as maps
import full_strand as fs
import schur_coordinates as sc
from shared_first import column_factor, SIZES
from coupled_first_replay import make_frame, PRIME as P, PLAN_SHA, RHS_SHA
from schur_side_images import cut

EDGE=(7,9,16,17)
SHARED='5a696be205573ddab96692bf890977ba19cd046d86582db9b3b76a98bc253bc3'

def loadz(f):return json.loads(gzip.decompress(f.read_bytes()))
def arr(f,digest):
    assert fs.fsha(f)==digest;return np.load(f,allow_pickle=False,mmap_mode='r')

def scaled(c,v):return fs.product(c,v) if len(c) and c.shape[1] and v.size else np.zeros((len(c),v.shape[1]),dtype=np.int64)

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate
    from native_su_values_v1 import evaluate as generic
    from g2_full_classical_screen_v1 import contents,sector
    from sp_full_profiles import profile
    from su3_stream_rank_v3 import random_word
    rt=a.root;out=a.out;out.mkdir(exist_ok=True,parents=True);plan=rt/'classical_order1.json';assert fs.fsha(plan)==PLAN_SHA;p=json.loads(plan.read_bytes())['incoming'];mu=tuple(p['mu']);ch=sector(mu,10);alpha=sorted({l for n in (10,11) for c in contents(mu,n) for l in c});profiles={i:maps.pattern(r) for i,r in enumerate(p['sources'])};ts={i:maps.pattern(t) for i,t in enumerate(p['targets'])};cs=defaultdict(list)
    for c in contents(mu,10):cs[profile(c)].append(c)
    assert set(EDGE)=={s for s in profiles if 10 in p['sources'][s]['targets']}
    f={i:make_frame(t,alpha,12,88304773+229*i) for i,t in enumerate(p['targets'])};cuts={};train=[];hold=[];rhs=[];nn=0
    for i,frame in f.items():
        cuts[i]=slice(nn,nn+frame.ncols);train.extend(range(nn,nn+8*frame.cpp));hold.extend(range(nn+8*frame.cpp,nn+frame.ncols));nn+=frame.ncols
        if i in (6,10,13,14,15,16):
            d=rt/('joint/edges' if i in (6,10) else 'first_profile16/coupled');meta=json.loads((d/f'fresh_{i}.json').read_bytes());b=arr(d/f'fresh_{i}.npy',meta['vectors_sha256']);b=b if i in (6,10) else b[2];assert meta['seed']==88304773+229*i and meta['prime']==P and (meta.get('formal_rhs_stream_sha256') or meta.get('rhs_stream_sha256'))==RHS_SHA and meta['residual_sha256']==fs.harray(b)
        else: assert str(i) not in p['first_rhs_trace_terms'] and not any(i in p['sources'][s]['targets'] for s in (19,20,21,22));b=np.zeros(frame.ncols,dtype='<u4')
        rhs.extend(b)
    rhs=fs.packed(rhs);assert nn==2244 and len(train)==1496 and len(hold)==748
    controls={}; hashes={}
    def images(s,words):
        dest=np.zeros((len(words),nn),dtype='<u4');support=p['sources'][s]['targets'];pat=profiles[s]
        for j in range(0,len(words),24):
            es=[project(alg.parity_project({w:1},1),ch,P) for w in words[j:j+24]]
            for t in support:
                v=evaluate(es,f[t],pat,ts[t]);dest[j:j+len(es),cuts[t]]=v
                if (s,t) not in controls:
                    for k,e in enumerate(es):
                        if not np.any(v[k,:f[t].cpp]):continue
                        g=cut(f[t],0,1);d={w:int(z)%P for w,z in alg.differential(e).items() if z%P and profile(tuple(l for tr in w for l in tr))==ts[t]};np.testing.assert_array_equal(v[k,:g.cpp],generic([d],g)[0]);controls[s,t]=fs.harray(v[k,:g.cpp]);break
        return dest
    if a.mode=='edges':
        rng=random.Random(202611190518);orig=loadz(rt/'target10_selected/recipes.json.gz');assert fs.fsha(rt/'target10_selected/recipes.json.gz')=='f426c15ae047a535ef7e8b5babf594ef12a68ff0333505f352cb74cc1c399118';orig=[(s,maps.decode(w)) for s,w in orig];assert len(orig)==704;seen=set(orig);rec=list(orig)
        for s in EDGE:
            added=0
            while added<a.additional:
                w=random_word(rng.choice(cs[profiles[s]]),rng,6)
                if w is not None and (s,w) not in seen:seen.add((s,w));rec.append((s,w));added+=1
        y=np.zeros((len(rec),nn),dtype='<u4')
        for s in EDGE:
            ix=[i for i,(k,w) in enumerate(rec) if k==s];y[ix]=images(s,[rec[i][1] for i in ix]);print('source',s,len(ix),flush=True)
        assert fs.harray(y[:704])==a.old_matrix_data_hash
        recipe=out/'expanded_edges.json.gz';recipe.write_bytes(gzip.compress(json.dumps(rec,separators=(',',':')).encode(),mtime=0));outrec=dict(recipe_sha256=fs.fsha(recipe),old_rows=704,additional_per_source=a.additional,first_704_entries_sha256=fs.harray(y[:704]),rows_by_source=dict(Counter(s for s,w in rec)))
    elif a.mode=='kernel':
        assert a.cluster is not None;cluster=a.cluster;coords=cluster/'su3_j34_reduced_20261119/replay';bases=cluster/'su3_j34_reduced_20261119/inputs';bundled=rt/'other_maps/task_5/minor.json.gz';_,words,_=sc.m.read_inputs(plan,bundled);assert len(words)==1112
        labels=[(s,j) for s,n in SIZES for j in range(n)];assert len(labels)==4125; pv=np.zeros((1112,nn),dtype='<u4');pos=0
        for ix,s in enumerate(sc.m.SOURCES):
            ix0=[i for i,(k,_) in enumerate(words) if k==ix];pv[ix0]=images(s,[words[i][1] for i in ix0]);pos+=len(ix0)
        assert pos==1112;v=np.zeros((4125,nn),dtype='<u4');pos=0
        for s,n in SIZES[:4]:
            idx,ww,k,cols,bs=sc.source_basis(argparse.Namespace(bases=bases),p,words,s);cmeta=json.loads((coords/f'coordinates_{s}.json').read_bytes());assert cmeta['source_basis_file_sha256']==bs and cmeta['source']==s and cmeta['pivot_minor_sha256']==sc.m.INPUT
            c=arr(coords/f'coefficients_{s}.npy',cmeta['coefficient_file_sha256']);assert c.shape==(n,1112) and len(ww[k:])==n
            if n:v[pos:pos+n]=(images(s,ww[k:]).astype(np.int64)-scaled(c,pv))%P
            pos+=n
        assert pos==3900;_,ww=fs.readwords(rt/'other_maps/task_7/minor.json.gz',fs.EXTRA);assert len(ww)==225;v[pos:]=images(6,ww);strand_record={}
        for t,(task,s,sha,_) in fs.STRANDS.items():
            d=cluster/f'su3_j34_full_strand_20261119/output/target_{t}';tr=json.loads((d/'result.json').read_bytes());assert tr['prime']==P and tr['minor_sha256']==sha and tr['strand_source']==s
            _,ww=fs.readwords(rt/f'strands/task_{task}/minor.json.gz',sha);c=arr(d/'coefficients.npy',tr['coefficient_file_sha256']);ns=len(c);assert c.shape==(4125 if t==9 else 3900,len(ww));vv=images(s,ww);delta=scaled(c,vv);v[:ns]=(v[:ns].astype(np.int64)-delta)%P;strand_record[str(s)]={str(u):int(np.count_nonzero(delta[:,cuts[u]])) for u in p['sources'][s]['targets']};print('strand',s,len(ww),strand_record[str(s)],flush=True)
        assert set(strand_record['18'])=={'11','12','14'};d=cluster/'su3_j34_shared_rows_20261119/zero';rf=d/'result.json';assert fs.fsha(rf)==SHARED;r=json.loads(rf.read_bytes());cc=arr(d/'shared_coefficients.npy',r['coefficients_sha256']);ii={k:i for i,k in enumerate(labels)};sel=[ii[tuple(z)] for z in r['pivot_source_rows']];free=[ii[tuple(z)] for z in r['free_source_rows']];assert sorted(sel+free)==list(range(4125)) and len(sel)==2606 and len(free)==1519;np.testing.assert_array_equal(cc[sel],np.eye(2606,dtype='<u4'))
        y=fs.packed((v[free].astype(np.int64)-scaled(cc[free],v[sel]))%P);chk={str(t):int(np.count_nonzero(y[:,cuts[t]])) for t in f};assert not chk['10'];outrec=dict(shared_recipe_receipt_sha256=SHARED,shared_coefficients_sha256=r['coefficients_sha256'],free_labels=[labels[i] for i in free],nonzero_on_fresh_coordinates=chk,strand_side_nonzeros=strand_record)
    else:
        y0=[];parents={}
        for k in ('edges','kernel'):
            rr=json.loads((out/f'{k}.json').read_bytes());assert rr['mode']==k and rr['plan_sha256']==PLAN_SHA and rr['rhs_sha256']==fs.harray(rhs);vv=arr(out/f'{k}.npy',rr['file_sha256']);y0.append(vv);parents[k]=rr['file_sha256']
        y=np.vstack(y0);rr,cc,dt,_=column_factor(y[:,train],panel=64);cols=[train[k] for k in cc];minor=nmod_mat(y[np.ix_(rr,cols)].tolist(),P);assert dt==int(minor.det()) and dt
        co=np.zeros(len(y),dtype='<u4');co[rr]=[int(z[0]) for z in minor.transpose().solve(nmod_mat([[int(rhs[i])] for i in cols],P)).tolist()];res=fs.packed((rhs.astype(np.int64)-scaled(co[None,:],y)[0])%P);assert not res[cols].any();bad=[j for j in train+hold if res[j]];ad=None
        if bad:
            j=bad[0];bcols=cols+[j];ad=int(nmod_mat(np.vstack([y[np.ix_(rr,bcols)],rhs[bcols]]).tolist(),P).det());assert ad==dt*int(res[j])%P and ad
        rec=dict(mode='merge',prime=P,plan_sha256=PLAN_SHA,source_matrix_hashes=parents,source_rows=len(y),training_rank=len(rr),training_width=len(train),withheld_width=len(hold),minor_determinant=dt,augmented_selected_minor=ad,first_extra_column=(bad or [None])[0],pivot_rows=rr,pivot_columns=cols,coefficient_count=int(np.count_nonzero(co)),training_nonzero=int(np.count_nonzero(res[train])),withheld_nonzero=int(np.count_nonzero(res[hold])),profiles={str(t):{'width':g.ncols,'training_nonzero':int(np.count_nonzero(res[cuts[t]][:8*g.cpp])),'withheld_nonzero':int(np.count_nonzero(res[cuts[t]][8*g.cpp:]))} for t,g in f.items()},rhs_sha256=fs.harray(rhs),coefficient_file_sha256=fs.write_array(out/'solution.npy',co),residual_file_sha256=fs.write_array(out/'residual.npy',res),residual_sha256=fs.harray(res));fs.putjson(out/'merge.json',rec);print('merge',len(y),len(rr),rec['training_nonzero'],rec['withheld_nonzero'],flush=True);return
    outrec.update(mode=a.mode,prime=P,plan_sha256=PLAN_SHA,frame_points=12,rows=len(y),width=nn,rhs_sha256=fs.harray(rhs),file_sha256=fs.write_array(out/f'{a.mode}.npy',y),differential_controls={f'{s}:{t}':v for (s,t),v in controls.items()});fs.putjson(out/f'{a.mode}.json',outrec);print(a.mode,len(y),outrec['file_sha256'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['edges','kernel','merge']);ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--root',type=Path,default=maps.HERE);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--cluster',type=Path);ap.add_argument('--additional',type=int,default=512);ap.add_argument('--old-matrix-data-hash',default='a2595f9a17b4700c030bc8704d6ef261007463e5a687068a65f4801b7d6dde3f');a=ap.parse_args();assert __debug__ and a.additional>0;run(a)
