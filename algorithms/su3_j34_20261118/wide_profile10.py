"""Evaluate an overdetermined projection and certify its row span over F_1000003.

The first 280 points are the archived profile-10 frame. Append 420 points at
seed 384116827. The 1,519 carried directions and the baseline are identically
zero on this profile by source incidence. No all-profile lifting claim follows.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import four_source_refit as refit

f, fs, P = refit.f, refit.fs, refit.P


def pivots(a):
    return f.joint.pivot_indices(f.matrix(a))


def certify(A, b, panel=64):
    """Schur elimination with dependent rows allowed in each panel; explicit dual/primal."""
    n,w=A.shape; assert b.shape==(w,) and n and w and min(A.min(),b.min())>=0 and max(A.max(),b.max())<P
    z=np.vstack((A,b)).astype('<u4'); U=[]; I=[]; J=[]; inverses=[]; blocks=[]; det=1
    for start in range(0,n,panel):
        stop=min(n,start+panel); orig=z[start:stop]; ii=pivots(orig.T)
        if not ii:assert not orig.any();continue
        xx=orig[ii]; jj=pivots(xx); assert len(ii)==len(jj) and set(J).isdisjoint(jj)
        sq=f.matrix(xx[:,jj]); value=int(sq.det()); assert value
        inv=fs.packed(sq.inv().tolist()); norm=fs.packed(fs.product(inv,xx))
        np.testing.assert_array_equal(norm[:,jj],np.eye(len(ii),dtype='<u4'))
        np.testing.assert_array_equal(fs.product(orig[:,jj],norm),orig)
        for k in range(stop,n+1,384):
            t=z[k:min(k+384,n+1)];t[:]=(t.astype(np.int64)-fs.product(t[:,jj],norm))%P
        U.append(norm);inverses.append(inv);I.extend(start+i for i in ii);J.extend(jj);blocks.append((jj,value));det=det*value%P
    rem=z[-1];bad=np.flatnonzero(rem);assert not rem[J].any() and det
    info=dict(rank=len(I),augmented_rank=len(I)+int(len(bad)>0),pivot_rows=I,pivot_columns=J,minor_determinant=det,panel_orders=[len(jj) for jj,d in blocks],panel_determinants=[d for jj,d in blocks],matrix_shape=[n,w],matrix_sha256=fs.harray(A),rhs_sha256=fs.harray(b))
    if len(bad):
        k=int(bad[0]); dual=np.zeros(w,dtype=np.int64);dual[k]=1
        for u,(jj,_) in zip(reversed(U),reversed(blocks)):
            assert not dual[jj].any();dual[jj]=(-u.astype(np.int64)@dual)%P
        for i in range(0,n,128):assert not (A[i:i+128].astype(np.int64)@dual%P).any()
        pair=int(b.astype(np.int64)@dual%P);assert pair==int(rem[k])!=0
        info.update(certificate='separator',extra_column=k,annihilated_active_rows=n,pairing=pair,augmented_minor_determinant=det*pair%P)
        return info,fs.packed(dual)
    # Solve the transposed system on the chosen square minor by processing
    # row panels backwards. Extract original coefficients by forward substitution
    # on transposed pivots (the inactive source rows receive coefficient zero).
    if I:
        M=f.matrix(A[np.ix_(I,J)]); cs=M.transpose().solve(f.matrix(b[J,None]));co=np.zeros(n,dtype=np.int64);co[I]=[int(row[0]) for row in cs.tolist()]
    else:co=np.zeros(n,dtype=np.int64)
    np.testing.assert_array_equal(refit.multiply(co,A),b);info['certificate']='coefficients';return info,fs.packed(co)


def main(a):
    previous=a.root/'out' if a.previous is None else a.previous
    fsargs=argparse.Namespace(root=a.root,previous=previous,runtime=a.runtime,workers=a.workers)
    sys.path.insert(0,str(a.runtime)); fit,built,prior,p,words,_=refit.frozen(fsargs)
    assert a.points==420 and 10*a.points+2800>len(words)==6848
    path=a.large/'matrix.npy'; active=fs.packed(np.load(path,allow_pickle=False));base=fs.packed(np.load(a.large/'rhs.npy',allow_pickle=False));assert fs.harray(active)==prior['matrix_data_sha256'] and fs.harray(base)==prior['rhs_data_sha256']
    path=previous/'new_rows.npy';assert fs.fsha(path)==built['new_matrix_sha256']; new=np.load(path,mmap_mode='r',allow_pickle=False); assertshape=(4096,5044);assert new.shape==assertshape
    mtrx=a.out/'matrix.npy';right=a.out/'rhs.npy'; assert not mtrx.exists() and not right.exists();a.out.mkdir(parents=True,exist_ok=True)
    H,y=refit.frame_matrix(fsargs,p,words,a.points,refit.WIDE_SEED);assert not H[2752:4271].any()
    A=np.concatenate((np.vstack((active,new[:,:2800])),np.vstack((H[:2752],H[4271:]))),axis=1);b=np.r_[base,y];del H,active,new
    assert A.shape==(6848,7000);fs.write_array(mtrx,A);fs.write_array(right,b);print('constructed',A.shape,flush=True)
    rec,cert=certify(A,b);sha=fs.write_array(a.out/(rec['certificate']+'.npy'),cert)
    rec.update(prime=P,source_profiles=list(f.SOURCES),effective_words=len(words),all_directions=8367,structurally_zero=1519,target_profile=10,target_dimension=p['targets'][10]['dimension'],old_points=280,old_seed=prior['seed'],old_width=2800,new_points=a.points,new_seed=refit.WIDE_SEED,new_width=len(y),prior_fit_matrix_sha256=fit['joined_rows_sha256'],new_words_sha256=built['new_words_sha256'],old_recipes_sha256=f.v.RECIPES_SHA,certificate_sha256=sha)
    fs.putjson(a.out/'wide.json',rec);print('ranks',rec['rank'],rec['augmented_rank'],rec['certificate'],rec.get('pairing'),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--previous',type=Path);ap.add_argument('--large',type=Path,required=True);ap.add_argument('--runtime',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--workers',type=int,default=32);ap.add_argument('--points',type=int,default=420);a=ap.parse_args();assert __debug__ and a.workers>0;main(a)
