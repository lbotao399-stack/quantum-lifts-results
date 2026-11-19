"""Extend a frozen rank certificate on the SAME 8,367 correction directions.

First replay the full original matrix product and the 120 withheld coordinates.
Then compute both ranks on the old kernel. A different seed is held out of that solve.
All ranks and coefficients are in F_1000003. Run with assertions enabled.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np
import four_source_joint as f

v, fs, P = f.v, f.fs, f.P
FRESH_SEED = 713904227
WIDE_SEED = 384116827
FROZEN = {'fit.json': '2eda7b644876fc683a51f5889ce928848de1cb262120814ed4737266cb9c5d1e',
          'validation.json': '081d7ad7c3bd574685b9dbb60ab79d50cc9ee5637cf1cf945a88c9486248b217',
          'build.json': 'b13601292b424c2a6b145a7fd632a7650d2b138350f17369f71e35ebd6301388',
          'coefficients.npy': '48594e1248354af2a887fdf031ab6072e1300b7ed81c318f231338a30aa3bbe4',
          'new_words.json.gz': 'b984cb977c1ae09ca97b65194b6cb3b037bb29900971a60a7373b94ebee5b87b'}


def frozen(a):
    assert all(fs.fsha(a.previous/k) == h for k,h in FROZEN.items())
    fit = json.loads((a.previous/'fit.json').read_bytes())
    build = json.loads((a.previous/'build.json').read_bytes())
    prior = json.loads((a.root/'allcoord/profile10_span.json').read_bytes())
    assert fs.fsha(a.root/'allcoord/profile10_span.json') == f.SPAN
    p, words, _ = v.inputs(a.root, a.root/'allcoord'); words += f.read_words(a.previous/'new_words.json.gz')
    assert len(words) == len(set(words)) == 6848
    c = np.load(a.previous/'coefficients.npy',allow_pickle=False).astype(np.int64)
    assert len(c)==8367 and c.min()>=0 and c.max()<P
    return fit,build,prior,p,words,c


def assembled(a,fit,build):
    A=np.zeros((8367,5044),dtype='<u4'); y=np.empty(5044,dtype='<u4')
    for name,part in [('edges',slice(0,2752)),('kernel',slice(2752,4271))]:
        path=a.archive/(name+'.npy'); assert fs.fsha(path)==f.joint.FILES[name]
        A[part,2800:]=np.load(path,allow_pickle=False)
    A[:2752,:2800]=np.load(a.large/'matrix.npy',allow_pickle=False)
    y[:2800]=np.load(a.large/'rhs.npy',allow_pickle=False); y[2800:]=f.joint.rhs(a.root)
    path=a.previous/'new_rows.npy'; assert fs.fsha(path)==build['new_matrix_sha256']
    A[4271:]=np.load(path,allow_pickle=False)
    assert fs.harray(A)==fit['joined_rows_sha256'] and fs.harray(y)==fit['joined_rhs_sha256']
    assert A.max()<P and y.max()<P
    return A,y


def multiply(c,A):
    # A different arithmetic path from the elimination: integer matrix-vector product.
    assert len(c)*(P-1)**2<2**63 and c.min()>=0 and c.max()<P
    y=np.zeros(A.shape[1],dtype=np.int64)
    for i in range(0,len(A),192):y=(y+c[i:i+192]@A[i:i+192].astype(np.int64))%P
    return y


def init(runtime,p,points,seed):
    global frames,charge,alg,project,derivative,evaluate,tpat,plan
    sys.path.insert(0,str(runtime)); import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from native_su_values_v1 import evaluate
    plan=p; fr,charge=v.frame(p,points,seed); tpat=v.maps.pattern(p['targets'][10]); frames=[f.span.part(fr,i,min(i+16,fr.P)) for i in range(0,fr.P,16)]


def block(task):
    kind,indices,rows=task
    if kind=='rhs':return kind,indices,fs.packed(np.concatenate([evaluate([dict(rows)],g)[0] for g in frames]))
    s=kind; ex=[project(alg.parity_project({w:1},1),charge,P) for w in rows]; pat=v.maps.pattern(plan['sources'][s]);return s,indices,fs.packed(np.concatenate([derivative(ex,g,pat,tpat) for g in frames],axis=1))


def frame_matrix(a,p,words,points,seed):
    def tasks():
        for batch in v.batches(list(v.rhs_terms(a.root,p).items()),192):yield 'rhs',None,batch
        yield from f.tasks_for(words,8)
    A=np.zeros((8367,10*points),dtype='<u4'); y=np.zeros(10*points,dtype=np.int64);done=np.zeros(8367,dtype=bool)
    args=(a.runtime,p,points,seed)
    with ProcessPoolExecutor(max_workers=a.workers,initializer=init,initargs=args) as pool:
        for kind,ids,value in pool.map(block,tasks()):
            if kind=='rhs':y=(y+value)%P;continue
            idx=[i if i<2752 else i+1519 for i in ids]
            assert not done[idx].any() and value.shape==(len(ids),len(y));done[idx]=True;A[idx]=value
    assert done[:2752].all() and not done[2752:4271].any() and done[4271:].all() and int(A.max())<P
    return A,fs.packed(y)


def solve_columns(T,B):
    # Factor just once. dot products in fs.product are split below 2**53.
    inv=fs.packed(T.inv().tolist());return fs.packed(fs.product(inv,B)),inv


def extension(A,y,H,u,c,I,J,det):
    # Previous rank r is verified by exact reconstruction on every discarded old column.
    assert not (multiply(c,A)-y).any()
    r=len(I);assert len(set(I))==len(set(J))==r
    used=set(J);free=[j for j in range(A.shape[1]) if j not in used]
    T=f.matrix(A[np.ix_(I,J)]); assert int(T.det())==det!=0
    legs=np.concatenate((A[np.ix_(I,free)],H[I]),axis=1);B,inv=solve_columns(T,legs)
    Z=np.empty(H.shape,dtype='<u4');nf=len(free)
    for i in range(0,len(A),160):
        sl=slice(i,min(i+160,len(A))); rows=A[sl][:,J]; rec=fs.product(rows,B)
        np.testing.assert_array_equal(rec[:,:nf],A[sl][:,free]); Z[sl]=(H[sl].astype(np.int64)-rec[:,nf:])%P
    assert not Z[I].any(); err=(u.astype(np.int64)-multiply(c,H))%P
    n=int(f.matrix(Z).rank());aug=int(f.matrix(np.vstack((Z,err))).rank());assert aug in (n,n+1)
    ii=f.joint.pivot_indices(f.matrix(Z).transpose()) if n else []; jj=f.joint.pivot_indices(f.matrix(Z[ii])) if n else []
    assert len(ii)==len(jj)==n
    minor=f.matrix(Z[np.ix_(ii,jj)]) if n else None;md=int(minor.det()) if n else 1; assert md
    h=(np.array([int(k[0]) for k in minor.transpose().solve(f.matrix(err[jj,None])).tolist()],dtype=np.int64) if n else np.array([],dtype=np.int64));rem=(err-h@Z[ii].astype(np.int64))%P
    assert bool(rem.any())==(aug>n)
    out=dict(rank=r+n,augmented_rank=r+aug,old_rank=r,old_free_columns_reproduced=nf,old_rows_reproduced=len(A),old_minor_determinant=int(det),extra_rank=n,extra_augmented_rank=aug,extra_minor_determinant=md,extra_pivot_rows=ii,extra_pivot_columns=jj,full_minor_determinant=int(det)*md%P,member=(n==aug))
    if n==aug:
        # e_i - A_i[J] T^-1 e_I annihilates the old frame.
        d=c.copy();d[ii]=(d[ii]+h)%P
        q=(h@A[np.ix_(ii,J)].astype(np.int64))%P;d[I]=(d[I]-fs.product(fs.packed(q[None,:]),inv)[0])%P
        np.testing.assert_array_equal(multiply(d,A),y);np.testing.assert_array_equal(multiply(d,H),u)
        out.update(old_residual_nonzero=0,new_residual_nonzero=0,coefficient_nonzeros=int(np.count_nonzero(d)))
        return out,'coefficients',fs.packed(d)
    k=int(np.flatnonzero(rem)[0]);t=np.zeros(len(err),dtype=np.int64);t[k]=1
    if n:t[jj]=[-int(z[0]) for z in minor.solve(f.matrix(Z[ii,k,None])).tolist()];t%=P
    dual=np.zeros(A.shape[1]+len(err),dtype=np.int64);dual[A.shape[1]:]=t;dual[J]=(-B[:,nf:].astype(np.int64)@t)%P
    for i in range(0,len(A),256):assert not ((A[i:i+256].astype(np.int64)@dual[:len(y)]+H[i:i+256].astype(np.int64)@t)%P).any()
    pair=int((y.astype(np.int64)@dual[:len(y)]+u.astype(np.int64)@t)%P);assert pair==int(rem[k])!=0
    out.update(annihilated_rows=len(A),pairing=pair,extra_column=k);return out,'separator',fs.packed(dual)


def main(a):
    fit,build,prior,p,words,c=frozen(a); A,y=assembled(a,fit,build);np.testing.assert_array_equal(multiply(c,A),y)
    diag=json.loads((a.previous/'validation.json').read_bytes()); assert diag['seed']==f.VALIDATION_SEED
    H,u=frame_matrix(a,p,words,12,diag['seed']);value=multiply(c,H);err=(u.astype(np.int64)-value)%P
    for k,x in [('rhs',u),('lhs',value),('residual',err)]:np.testing.assert_array_equal(x,diag[k]);
    assert fs.harray(err)==diag['residual_sha256'] and np.count_nonzero(err)==120
    print('replayed 5044 fitted and 120 withheld coordinates',flush=True)
    stage=fit['stages']['extension'];I=prior['pivot_rows']+stage['quotient_pivot_rows'];J=prior['pivot_columns']+stage['quotient_pivot_columns']; assert len(I)==len(J)==stage['rank']==4108
    rec,kind,d=extension(A,y,H,u,c,I,J,prior['pivot_minor_determinant']*stage['quotient_minor_determinant']%P);path=a.out/(kind+'.npy');sha=fs.write_array(path,d)
    rec.update(prime=P,rows=len(A),fitted_width=len(y),appended_width=len(u),original_rhs_sha256=fs.harray(y),original_matrix_sha256=fs.harray(A),original_fitted_residual_nonzero=0,archived_solution_sha256=FROZEN['coefficients.npy'],appended_matrix_sha256=fs.harray(H),appended_rhs_sha256=fs.harray(u),original_withheld_residual_nonzero=120,original_withheld_residual_sha256=fs.harray(err),certificate_type=kind,certificate_file_sha256=sha,precommitted_fresh_seed=FRESH_SEED)
    fs.putjson(a.out/'rank.json',rec);print('refit ranks',rec['rank'],rec['augmented_rank'],flush=True)
    if kind=='coefficients':
        assert FRESH_SEED!=diag['seed'];K,w=frame_matrix(a,p,words,12,FRESH_SEED);assert fs.fsha(path)==sha;d=np.load(path,allow_pickle=False).astype(np.int64);val=multiply(d,K);res=(w.astype(np.int64)-val)%P
        fs.putjson(a.out/'fresh.json',dict(seed=FRESH_SEED,prime=P,profile=10,points=12,width=len(w),rhs=w.tolist(),lhs=val.tolist(),residual=res.tolist(),nonzero=int(np.count_nonzero(res)),residual_sha256=fs.harray(res),frozen_coefficients_sha256=sha));print('new holdout',np.count_nonzero(res),'of',len(res),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('root','previous','archive','large','out','runtime'):ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--workers',type=int,default=16);a=ap.parse_args();assert __debug__ and a.workers>0;a.out.mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(a.runtime.resolve()));main(a)
