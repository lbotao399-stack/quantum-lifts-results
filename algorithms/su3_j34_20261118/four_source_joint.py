"""Exact four-source extension of the two archived SU(3) j=34 evaluation frames.

Each projected trace has ONE coefficient on all of its target profiles. All 1519
carried rows are included. The held-out seed is evaluated only after the fit.
"""
import argparse
import gzip
import json
import sys
import random
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from flint import nmod_mat
import all_coordinate_rank as joint
import full_strand as fs
import profile10_span as span
import profile10_validation as v

P = v.PRIME
SOURCES = v.EDGE
ORDER = {7: [2,6,10], 9: [3,6,8,10,11], 16: [10,13], 17: [10,11,13,14]}
SPAN = 'ae565ce3f91869bac4dfbb026dda9b185e23591490af4995970bbecd557814d1'
WORD_SEED = 790231487
VALIDATION_SEED = 925173601


def write_words(path, words):
    x = json.dumps(words, separators=(',',':')).encode()
    path.write_bytes(gzip.compress(x, mtime=0)); return fs.fsha(path)


def read_words(path):
    return [(s,v.maps.decode(w)) for s,w in json.loads(gzip.decompress(path.read_bytes()))]


def make_words(root,runtime,per):
    sys.path.insert(0,str(runtime))
    from g2_full_classical_screen_v1 import contents
    from sp_full_profiles import profile
    from su3_stream_rank_v3 import random_word
    p,archived,_=v.inputs(root,root/'allcoord'); assert {s:p['sources'][s]['targets'] for s in SOURCES}==ORDER
    data=json.loads((root/'allcoord/profile10_external.json').read_bytes()); assert data['certificate_sha256']==SPAN
    start=[(c['source_profile'],v.maps.decode(c['recipe'])) for c in data['cases']]; assert [s for s,w in start]==list(SOURCES)
    seen=set(archived); assert not seen.intersection(start); seen.update(start)
    cs=contents(tuple(p['mu']),p['source_length']); choices={s:[w for w in cs if profile(w)==v.maps.pattern(p['sources'][s])] for s in SOURCES}
    out=list(start); r=random.Random(WORD_SEED)
    for _ in range(1,per):
        for s in SOURCES:
            while True:
                w=random_word(r.choice(choices[s]),r,8)
                if w is not None and (s,w) not in seen:break
            out.append((s,w));seen.add((s,w))
    assert dict(Counter(s for s,w in out))=={s:per for s in SOURCES}
    return out


def initialize(root,rt,big_points=280):
    global frames,cuts,width,prime_ref
    sys.path.insert(0,str(rt)); from g2_full_classical_screen_v1 import contents
    prime_ref=json.loads((root/'allcoord/profile10_span.json').read_bytes()); assert fs.fsha(root/'allcoord/profile10_span.json')==SPAN
    span.setup(rt,root,root/'allcoord',big_points,prime_ref['seed']); p=span.plan
    alpha=sorted({l for n in (p['source_length'],p['target_length']) for c in contents(tuple(p['mu']),n) for l in c})
    frames={t:v.make_frame(rec,alpha,12,88304773+229*t) for t,rec in enumerate(p['targets'])}
    start=span.fr.ncols;cuts={}
    for t,f in frames.items():cuts[t]=slice(start,start+f.ncols);start+=f.ncols
    width=start;assert width==span.fr.ncols+2244


def ev(task):
    s,ix,ww=task; p=span.plan;pat=v.maps.pattern(p['sources'][s]); assert p['sources'][s]['targets']==ORDER[s]
    es=[span.project(span.alg.parity_project({w:1},1),span.charge,P) for w in ww];A=np.zeros((len(ww),width),dtype='<u4')
    A[:,:span.fr.ncols]=np.concatenate([span.derivative(es,f,pat,span.tpat) for f in span.frames],axis=1)
    for t in ORDER[s]:A[:,cuts[t]]=span.derivative(es,frames[t],pat,v.maps.pattern(p['targets'][t]))
    return ix,A


def tasks_for(words,n=6):
    for s in SOURCES:
        idx=[i for i,(k,w) in enumerate(words) if k==s]
        for j in range(0,len(idx),n):
            ids=idx[j:j+n];yield s,ids,[words[i][1] for i in ids]


def check_rows(root,A,words,npoints):
    cases=json.loads((root/'allcoord/profile10_external.json').read_bytes())['cases']; checks={}
    from sp_full_profiles import profile
    for i,row in enumerate(cases):
        s,w=words[i]; assert s==row['source_profile'] and w==v.maps.decode(row['recipe'])
        if npoints==280:
            column=np.asarray(prime_ref['annihilator'],dtype=np.int64)
            assert fs.harray(A[i,:2800])==row['evaluations_sha256'] and int(A[i,:2800]@column%P)==row['pairing']
        exp=span.project(span.alg.parity_project({w:1},1),span.charge,P); dx=span.alg.differential(exp)
        for t in ORDER[s]:
            target=v.maps.pattern(span.plan['targets'][t]); terms={u:int(k)%P for u,k in dx.items() if k%P and profile(tuple(x for tr in u for x in tr))==target}
            g=span.part(frames[t],0,1); out=np.zeros(g.ncols,dtype=np.int64); items=list(terms)
            for u in v.batches(items,24):out=(out+np.array([terms[h] for h in u],dtype=np.int64)@g.eval_words(u,dtype=np.int64))%P
            np.testing.assert_array_equal(out,A[i,cuts[t]][:g.ncols]);checks[f'{s}:{t}']=fs.harray(out)
    return checks


def build(args):
    d=args.out;root=args.root;rt=args.runtime;words=make_words(root,rt,args.per_source); wh=write_words(d/'new_words.json.gz',words)
    initialize(root,rt,args.big_points);A=np.empty((len(words),width),dtype='<u4');done=np.zeros(len(words),dtype=bool)
    if args.workers==1: results=map(ev,tasks_for(words)); pool=None
    else:
        pool=ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(root,rt,args.big_points)); results=pool.map(ev,tasks_for(words))
    try:
        for i,(ix,block) in enumerate(results):
            assert not done[ix].any() and block.shape==(len(ix),width);A[ix]=block;done[ix]=True
            if (i+1)%100==0:print('new word blocks',i+1,flush=True)
    finally:
        if pool is not None:pool.shutdown()
    assert done.all() and A.max()<P; checks=check_rows(root,A,words,args.big_points)
    rec=dict(word_seed=WORD_SEED,per_source=args.per_source,new_words_sha256=wh,new_rows=len(words),large_points=args.big_points,shape=list(A.shape),source_dimensions={s:span.plan['sources'][s]['dimension'] for s in SOURCES},incidence=ORDER,explicit_checks=checks,new_matrix_sha256=fs.write_array(d/'new_rows.npy',A))
    fs.putjson(d/'build.json',rec);return rec


def matrix(row):return nmod_mat(row.tolist(),P)


def inherited_quotient(A,b,prior):
    """Eliminate the known square minor once; keep original row coefficients."""
    I=prior['pivot_rows'];J=prior['pivot_columns'];assert len(I)==len(J)==prior['rank']
    B=A[I]; S=matrix(B[:,J]);det=int(S.det());assert det==prior['pivot_minor_determinant'];inv=fs.packed(S.inv().tolist());C0=set(J);other=np.array([i for i in range(A.shape[1]) if i not in C0]);D=B[:,other]
    rows=np.vstack((A,b));C=np.empty((len(rows),len(I)),dtype='<u4');Y=np.empty((len(rows),len(other)),dtype='<u4')
    for off in range(0,len(rows),128):
        sl=slice(off,min(off+128,len(rows)));c=fs.packed(fs.product(rows[sl][:,J],inv));C[sl]=c;Y[sl]=(rows[sl][:,other].astype(np.int64)-fs.product(c,D))%P
    assert not Y[I].any() and np.array_equal(C[I],np.eye(len(I),dtype='<u4'))
    return I,J,other,C,Y,inv,D


def evaluate_stage(A,b,I,J,other,C,Y,inv,D,n,detail):
    """Eliminate on the exact quotient; restore a primal or dual on the full frame."""
    Z=Y[:n];y=Y[-1];rank=int(matrix(Z).rank());aug=int(matrix(np.vstack((Z,y))).rank());result=dict(rows=n,rank=len(I)+rank,augmented_rank=len(I)+aug,member=rank==aug)
    assert aug in (rank,rank+1)
    if not detail:return result,None
    ii=joint.pivot_indices(matrix(Z).transpose()); zz=Z[ii]; jj=joint.pivot_indices(matrix(zz));assert rank==len(ii)==len(jj);T=matrix(zz[:,jj]);det=int(T.det());assert det
    v0=np.array([int(a[0]) for a in T.transpose().solve(matrix(y[jj,None])).tolist()],dtype=np.int64);c=np.zeros(n,dtype=np.int64);c[ii]=v0;c[I]=(c[I]+C[-1].astype(np.int64)-v0@C[ii].astype(np.int64))%P
    rebuild=np.zeros(len(b),dtype=np.int64)
    for off in range(0,n,256):rebuild=(rebuild+c[off:off+256]@A[off:off+256].astype(np.int64))%P
    res=(b.astype(np.int64)-rebuild)%P;assert not res[J].any() and not res[other[jj]].any() and bool(res.any())==(rank!=aug)
    result.update(quotient_minor_determinant=det,quotient_pivot_rows=ii,quotient_pivot_columns=[int(other[k]) for k in jj],trial_residual_sha256=fs.harray(res))
    if rank==aug:
        result.update(coefficient_nonzeros=int(np.count_nonzero(c)),coefficient_data_sha256=fs.harray(c),nonzero_by_row_block=[int(np.count_nonzero(a)) for a in (c[:2752],c[2752:4271],c[4271:])]);return result,fs.packed(c)
    j=int(np.flatnonzero(res[other])[0]);t=np.zeros(len(other),dtype=np.int64);t[j]=1;t[jj]=[-int(z[0]) for z in T.solve(matrix(zz[:,j,None])).tolist()];t%=P
    ext=np.zeros(len(b),dtype=np.int64);ext[other]=t;ext[J]=(-fs.product(inv,fs.packed((D.astype(np.int64)@t%P)[:,None]))[:,0])%P
    for k in range(0,n,192):assert not (A[k:min(n,k+192)].astype(np.int64)@ext%P).any()
    pair=int(b.astype(np.int64)@ext%P);assert pair==int(res[other[j]])!=0;result.update(annihilated_rows=n,pairing=pair,dual_data_sha256=fs.harray(ext));return result,fs.packed(ext)


def solve(args):
    d=args.out; root=args.root; assert args.big_points==280
    old=args.archive; big=args.large; rec=json.loads((d/'build.json').read_bytes());assert rec['large_points']==280 and fs.fsha(d/'new_rows.npy')==rec['new_matrix_sha256']
    v.inputs(root,root/'allcoord');prior=json.loads((root/'allcoord/profile10_span.json').read_bytes()); assert fs.fsha(root/'allcoord/profile10_span.json')==SPAN
    for key in ('edges','kernel'): assert fs.fsha(old/f'{key}.npy')==joint.FILES[key]
    e,k=(np.load(old/f'{key}.npy',allow_pickle=False) for key in ('edges','kernel'));x=np.load(big/'matrix.npy',allow_pickle=False);y=np.load(big/'rhs.npy',allow_pickle=False);y0=joint.rhs(root);a=np.load(d/'new_rows.npy',allow_pickle=False)
    assert fs.harray(x)==prior['matrix_data_sha256'] and fs.harray(y)==prior['rhs_data_sha256'] and (e.shape,k.shape,x.shape,y.shape)==((2752,2244),(1519,2244),(2752,2800),(2800,)); assert a.shape[1]==5044
    b=np.concatenate((y,y0));A=np.zeros((4271+len(a),len(b)),dtype='<u4');A[:2752,:2800]=x;A[:2752,2800:]=e;A[2752:4271,2800:]=k;A[4271:]=a;del e,k,x,a
    I,J,other,C,Y,inv,D=inherited_quotient(A,b,prior);print('inherited exact minor',len(I),flush=True)
    stages={}; solution=None
    for label,n,detail in [('archived',4271,False),('four_witnesses',4275,False),('extension',len(A),True)]:
        stages[label],solution=evaluate_stage(A,b,I,J,other,C,Y,inv,D,n,detail)
        large_only=Y[:n,other<2800];bb=Y[-1,other<2800]; assert large_only.shape[1]==109
        stages[label]['large_profile10_rank']=len(I)+int(matrix(large_only).rank())
        stages[label]['large_profile10_augmented_rank']=len(I)+int(matrix(np.vstack((large_only,bb))).rank())
        print(label,stages[label]['rank'],stages[label]['augmented_rank'],flush=True)
    typ='coefficients' if stages['extension']['member'] else 'separator'; cert=fs.write_array(d/(typ+'.npy'),solution)
    audit=dict(prime=P,width=len(b),large_profile10_width=2800,original_all_profile_width=2244,joined_rows_sha256=fs.harray(A),joined_rhs_sha256=fs.harray(b),inherited_minor_sha256=SPAN,inherited_minor_order=len(I),inherited_minor_determinant=prior['pivot_minor_determinant'],new_words_sha256=rec['new_words_sha256'],old_joint_files=joint.FILES,stages=stages,certificate_type=typ,certificate_file_sha256=cert,precommitted_validation_seed=VALIDATION_SEED)
    fs.putjson(d/'fit.json',audit);return audit


def validate(args):
    d=args.out; root=args.root;sys.path.insert(0,str(args.runtime));from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from native_su_values_v1 import evaluate
    import su_trace_adapter_v1 as alg
    fit=json.loads((d/'fit.json').read_bytes());assert fit['precommitted_validation_seed']==VALIDATION_SEED
    if fit['certificate_type']!='coefficients':return
    fl=d/'coefficients.npy';assert fs.fsha(fl)==fit['certificate_file_sha256'];c=np.load(fl,allow_pickle=False).astype(np.int64);p,old,_=v.inputs(root,root/'allcoord');ww=read_words(d/'new_words.json.gz');assert fs.fsha(d/'new_words.json.gz')==fit['new_words_sha256']; assert len(c)==4271+len(ww)
    cc=np.concatenate((c[:2752],c[4271:]));words=old+ww;f,ch=v.frame(p,args.validation_points,VALIDATION_SEED);pat=v.maps.pattern(p['targets'][10]);rhs=np.zeros(f.ncols,dtype=np.int64);lhs=rhs.copy()
    for batch in v.batches(list(v.rhs_terms(root,p).items()),192):rhs=(rhs+evaluate([dict(batch)],f)[0])%P
    for s,ix,seq in tasks_for(words,8):
        if not cc[ix].any():continue
        ex=[project(alg.parity_project({w:1},1),ch,P) for w in seq];lhs=(lhs+cc[ix]@derivative(ex,f,v.maps.pattern(p['sources'][s]),pat))%P
    err=(rhs-lhs)%P; assert fit['certificate_type']=='coefficients';rec=dict(prime=P,seed=VALIDATION_SEED,profile=10,points=args.validation_points,structurally_zero_archived_kernel=1519,structurally_zero_baseline=True,frozen_solution_sha256=fit['certificate_file_sha256'],rhs=rhs.tolist(),lhs=lhs.tolist(),residual=err.tolist(),nonzero=int(np.count_nonzero(err)),width=len(err),residual_sha256=fs.harray(err));fs.putjson(d/'validation.json',rec);print('frozen fit: fresh profile 10',rec['nonzero'],'of',len(err),flush=True)


if __name__=='__main__':
    pa=argparse.ArgumentParser(description=__doc__);pa.add_argument('mode',choices=['all','build','solve','validate']);pa.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);pa.add_argument('--runtime',type=Path,default=v.maps.RUNTIME);pa.add_argument('--archive',type=Path);pa.add_argument('--large',type=Path);pa.add_argument('--out',type=Path,required=True);pa.add_argument('--per-source',type=int,default=1024);pa.add_argument('--workers',type=int,default=16);pa.add_argument('--big-points',type=int,default=280);pa.add_argument('--validation-points',type=int,default=12);args=pa.parse_args();assert __debug__ and min(args.per_source,args.workers,args.big_points,args.validation_points)>0;args.runtime=args.runtime.resolve();sys.path.insert(0,str(args.runtime));args.out.mkdir(parents=True,exist_ok=True)
    if args.mode in ('build','all'):build(args)
    if args.mode in ('solve','all'):solve(args)
    if args.mode in ('validate','all'):validate(args)
