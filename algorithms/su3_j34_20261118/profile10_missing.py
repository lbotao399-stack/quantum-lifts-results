"""Extend the frozen 7000-coordinate image using its 365-dimensional cokernel.

The audit uses the saved pivot order but refactors its minor independently. A
full-rank quotient yields one explicit old+new recipe solution on this frame.
It is then tested at a seed fixed before elimination. No sampled rank upper
bounds the full physical image; coefficients are over the indicated prime.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import hashlib
import json
import random
import sys
import numpy as np
import four_source_joint as f
import four_source_refit as refit
import replay_wide_pairings as independent

fs,v,P=f.fs,f.v,f.P
WORD_SEED=293671843
FRESH_SEED=635729813
FRAME=((280,419870023),(420,384116827))


def signature(expr):
    items=sorted((w,int(c)%P) for w,c in expr.items() if int(c)%P)
    if not items:return None
    inv=pow(items[0][1],-1,P)
    return hashlib.sha256(json.dumps([(w,(c*inv)%P) for w,c in items],separators=(',',':')).encode()).digest()


def generate(root,previous,runtime,per):
    sys.path.insert(0,str(runtime));import su_trace_adapter_v1 as alg
    from g2_full_classical_screen_v1 import contents, sector
    from native_su_projection_v1 import project
    from sp_full_profiles import profile
    from su3_stream_rank_v3 import random_word
    p,old,_=v.inputs(root,root/'allcoord'); assert fs.fsha(previous/'new_words.json.gz')==refit.FROZEN['new_words.json.gz']
    old+=f.read_words(previous/'new_words.json.gz');used=set(old);assert len(used)==len(old)==6848
    ct=contents(tuple(p['mu']),p['source_length']);charge=sector(tuple(p['mu']),p['source_length'])
    choices={s:[x for x in ct if profile(x)==v.maps.pattern(p['sources'][s])] for s in f.SOURCES}
    def key(w):return signature(project(alg.parity_project({w:1},1),charge,P))
    seen={(s,key(w)) for s,w in old};old_n=len(seen);rng=random.Random(WORD_SEED);new=[];attempts=Counter()
    for _ in range(per):
        for s in f.SOURCES:
            for _try in range(10000):
                attempts[s]+=1;w=random_word(rng.choice(choices[s]),rng,8)
                if not w or (s,w) in used:continue
                k=key(w)
                if k is None or (s,k) in seen:continue
                used.add((s,w));seen.add((s,k));new.append((s,w));break
            else:raise RuntimeError('could not find a distinct projected trace')
    assert len(new)==len(set(new))==len(seen)-old_n==4*per
    return p,old,new,dict(attempts=dict(attempts),archived_recipe_count=len(old),archived_projective_keys_including_zero=old_n)


def rows(words,p,runtime,workers,points,seed):
    out=np.empty((len(words),points*10),dtype='<u4');done=np.zeros(len(out),dtype=bool)
    with ProcessPoolExecutor(max_workers=workers,initializer=refit.init,initargs=(runtime,p,points,seed)) as pool:
        for s,ix,a in pool.map(refit.block,f.tasks_for(words,8)):
            assert not done[ix].any() and a.shape==(len(ix),out.shape[1]);out[ix]=a;done[ix]=True
    assert done.all();return out


def kernel(old,rec):
    I,J=rec['pivot_rows'],rec['pivot_columns'];chosen=set(J);free=[j for j in range(old.shape[1]) if j not in chosen];r=len(I); assert r==len(J) and len(free)>0
    u=fs.packed(old[np.ix_(I,J+free)]);lower=np.zeros((r,r),dtype='<u4');panels=[];start=0;det=1
    for size,expected in zip(rec['panel_orders'],rec['panel_determinants']):
        stop=start+size;sl=slice(start,stop);T=f.matrix(u[sl,sl]);d=int(T.det());assert d==expected!=0
        det=det*d%P;inv=fs.packed(T.inv().tolist());u[sl,start:]=fs.product(inv,u[sl,start:]);assert np.array_equal(u[sl,sl],np.eye(size,dtype='<u4'))
        for i in range(stop,r,320):
            block=u[i:i+320,start:];c=block[:,:size].copy();lower[i:i+len(c),sl]=c
            block[:]=(block.astype(np.int64)-fs.product(c,u[sl,start:]))%P
        panels.append((start,stop,inv));start=stop
    assert start==r and det==rec['minor_determinant']
    right=u[:,r:].copy()
    for start,stop,inv in reversed(panels):
        if stop<r:right[start:stop]=(right[start:stop].astype(np.int64)-fs.product(u[start:stop,stop:r],right[stop:]))%P
    return I,J,free,right,u[:,:r],lower,panels


def reduce(x,J,free,right):
    return fs.packed((x[:,free].astype(np.int64)-fs.product(x[:,J],right))%P)


def coefficients(y,upper,lower,panels):
    a=y[None,:].astype(np.int64).copy()
    for i,j,_ in panels:
        if j<a.shape[1]:a[:,j:]=(a[:,j:]-fs.product(a[:,i:j],upper[i:j,j:]))%P
    for i,j,inv in reversed(panels):
        a[:,i:j]=fs.product(a[:,i:j],inv)
        if i:a[:,:i]=(a[:,:i]-fs.product(a[:,i:j],lower[i:j,:i]))%P
    return a[0]


def fresh(root,p,words,c,runtime,workers):
    chosen=[i for i,z in enumerate(c) if z]; seq=[words[i] for i in chosen];co=c[chosen].astype(np.int64); width=120
    H=rows(seq,p,runtime,workers,12,FRESH_SEED);lhs=refit.multiply(co,H);fr,_=v.frame(p,12,FRESH_SEED)
    from native_su_values_v1 import evaluate
    rhs=np.zeros(width,dtype=np.int64)
    for g in v.batches(list(v.rhs_terms(root,p).items()),192):rhs=(rhs+evaluate([dict(g)],fr)[0])%P
    res=(rhs-lhs)%P;return dict(seed=FRESH_SEED,width=width,nonzero=int(np.count_nonzero(res)),lhs=lhs.tolist(),rhs=rhs.tolist(),residual=res.tolist(),residual_sha256=fs.harray(res))


def main(a):
    sys.path.insert(0,str(a.runtime));dest=a.out;dest.mkdir(parents=True,exist_ok=True);assert not (dest/'result.json').exists()
    audit=independent.run(a.wide,a.wide);rec=json.loads((a.wide/'wide.json').read_bytes());p,old,new,wr=generate(a.root,a.previous,a.runtime,a.per);word_hash=f.write_words(dest/'words.json.gz',new)
    oldA=np.load(a.wide/'matrix.npy',mmap_mode='r',allow_pickle=False);b=np.load(a.wide/'rhs.npy',allow_pickle=False);dual=np.load(a.wide/'separator.npy',allow_pickle=False).astype(np.int64)
    H=fs.packed(np.concatenate([rows(new,p,a.runtime,a.workers,n,seed) for n,seed in FRAME],axis=1));fs.write_array(dest/'new_rows.npy',H);print('distinct projected recipes evaluated',len(H),flush=True)
    vals=H.astype(np.int64)@dual%P;I,J,F,K,U,L,panels=kernel(oldA,rec);assert len(F)==365;old_nonpivot=sorted(set(range(len(oldA)))-set(I));assert len(old_nonpivot)==213
    for off in range(0,len(oldA),192):assert not reduce(oldA[off:off+192],J,F,K).any()
    # The first quotient column is the earlier dual functional, with unit free coefficient.
    jj=F.index(rec['extra_column']);assert dual[F[jj]]==1 and np.count_nonzero(dual[F])==1;np.testing.assert_array_equal((-K[:,jj].astype(np.int64))%P,dual[J]);Y=reduce(H,J,F,K);y=reduce(b[None,:],J,F,K)[0];np.testing.assert_array_equal(vals,Y[:,jj]);assert int(y[jj])==audit['rhs_pairing'];fs.write_array(dest/'quotient.npy',Y);fs.write_array(dest/'quotient_rhs.npy',y)
    per={}; witnesses={}
    for s in f.SOURCES:
        idx=[i for i,(src,w) in enumerate(new) if src==s]; hit=[i for i in idx if vals[i]];per[s]={'sample':len(idx),'rank':int(f.matrix(Y[idx]).rank()),'nonzero_on_old_separator':len(hit)}
        if hit:
            i=hit[0]; witnesses[s]=dict(index=i,pairing=int(vals[i]),recipe=new[i][1],row_data_sha256=fs.harray(H[i]));fs.write_array(dest/('witness_'+str(s)+'.npy'),H[i])
    ii=f.joint.pivot_indices(f.matrix(Y).transpose());Z=Y[ii];cols=f.joint.pivot_indices(f.matrix(Z));T=f.matrix(Z[:,cols]);det=int(T.det());assert det!=0 and len(cols)==len(ii);above=int(f.matrix(np.vstack((Y,y))).rank());assert above in (len(ii),len(ii)+1)
    result=dict(prime=P,frames=FRAME,source_profiles=list(f.SOURCES),word_seed=WORD_SEED,precommitted_holdout_seed=FRESH_SEED,word_sampling=wr,word_file_sha256=word_hash,old_rank=len(I),old_determinant=rec['minor_determinant'],width=len(b),old_rows_reduced_to_zero=len(oldA),new_rows=len(H),source_stats=per,witnesses=witnesses,quotient_rank=len(ii),quotient_augmented_rank=above,quotient_determinant=det,quotient_pivot_rows=ii,quotient_pivot_cols=cols,free_coordinates=F,new_rank=len(I)+len(ii),augmented_rank=len(I)+above,new_minor_determinant=int(det)*rec['minor_determinant']%P,quotient_data_sha256=fs.harray(Y),quotient_rhs_data_sha256=fs.harray(y),new_rows_data_sha256=fs.harray(H),old_record_sha256=independent.RECORD)
    if above==len(ii):
        nc=np.zeros(len(H),dtype=np.int64);nc[ii]=[int(t[0]) for t in T.transpose().solve(f.matrix(y[cols,None])).tolist()];np.testing.assert_array_equal(refit.multiply(nc,Y),y);diff=(b.astype(np.int64)-refit.multiply(nc,H))%P;oldc=np.zeros(len(oldA),dtype=np.int64);oldc[I]=coefficients(diff[J],U,L,panels);np.testing.assert_array_equal(refit.multiply(oldc,oldA),diff);c=fs.packed(np.r_[oldc,nc]);result.update(certificate='coefficients',nonzero_old=int(np.count_nonzero(oldc)),nonzero_new=int(np.count_nonzero(nc)),coefficients_file_sha256=fs.write_array(dest/'coefficients.npy',c));fs.putjson(dest/'result.json',result);print('ranks',result['new_rank'],result['augmented_rank'],'checking fresh frame',flush=True);h=fresh(a.root,p,old+new,c,a.runtime,a.workers);h['frozen_coefficients_sha256']=result['coefficients_file_sha256'];fs.putjson(dest/'fresh.json',h);print('fresh nonzero',h['nonzero'],flush=True)
    else:
        v0=fs.packed(T.transpose().solve(f.matrix(y[cols,None])).transpose().tolist());e=(y.astype(np.int64)-fs.product(v0,Z)[0])%P;idx=int(np.flatnonzero(e)[0]);t=np.zeros(len(F),dtype=np.int64);t[idx]=1;t[cols]=[-int(x[0])%P for x in T.solve(f.matrix(Z[:,idx,None])).tolist()];d=np.zeros(len(b),dtype=np.int64);d[F]=t;d[J]=(-K.astype(np.int64)@t)%P
        for B in (oldA,H):
            for off in range(0,len(B),64):assert not (B[off:off+64].astype(np.int64)@d%P).any()
        pair=int(b.astype(np.int64)@d%P);assert pair==e[idx]!=0;result.update(certificate='separator',pairing=pair,separator_sha256=fs.write_array(dest/'separator.npy',d));fs.putjson(dest/'result.json',result)
    print('complete',result['new_rank'],result['augmented_rank'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ('root','previous','wide','runtime','out'):ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--workers',type=int,default=16);ap.add_argument('--per',type=int,default=192);args=ap.parse_args();assert __debug__ and args.per>0;main(args)
