"""Replay the five-frame certificates over F_1000003 without the producer's solver.

With --wide and --quotient pass the original two job output directories and use
--folder for the five-frame job output; these inputs are retained on A8. With
neither option only the compact arrays in the repository are needed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

P = 1000003
RECORD = 'f0cad202c283e4f1ef46833f6c2d5b2abaf028ddf3844e827685901db4725aee'
GROUPS = [('five_profiles','side_fit',0,312),('fresh10_only','fresh10_only',312,1112),('all_five','all_five',0,1112)]

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def hv(x):return hashlib.sha256(np.asarray(x,dtype='<u4').tobytes()).hexdigest()

def array(p,shape):
    a=np.load(p,mmap_mode='r',allow_pickle=False);assert a.shape==shape and a.dtype.str=='<u4' and int(a.max())<P
    return a

def eliminate(a):
    """Exact row echelon form; no producer pivots for the rank calculation."""
    a=np.array(a,dtype=np.int64);r=0; d=1
    for j in range(a.shape[1]):
        hits=np.flatnonzero(a[r:,j]);
        if not hits.size:continue
        k=r+int(hits[0]);
        if k!=r:a[[k,r]]=a[[r,k]];d=-d
        v=int(a[r,j]);d=d*v%P
        t=a[r+1:,j]*pow(v,-1,P)%P;nz=np.flatnonzero(t)
        if len(nz):a[r+1+nz,j:]=(a[r+1+nz,j:]-t[nz,None]*a[r,j:])%P
        r+=1
        if r==len(a):break
    return r,d%P

def read(folder):
    assert __debug__ and sha(folder/'result.json')==RECORD
    rec=json.loads((folder/'result.json').read_bytes());assert (rec['prime'],rec['recipes'])==(P,7616)
    cert={}
    for key,name,lo,hi in GROUPS:
        p=folder/(name+'.npy');assert sha(p)==rec[key]['certificate_file_sha256'];cert[name]=array(p,(7616 if lo==0 and hi==312 else 7000+hi-lo,)).astype(np.int64)
    return rec,cert

def compact(folder):
    rec,cert=read(folder);z=array(folder/'schur.npy',(616,1112));u=array(folder/'schur_rhs.npy',(1112,))
    assert (hv(z),hv(u))==(rec['schur_sha256'],rec['schur_rhs_sha256']); report={}
    for key,name,lo,hi in GROUPS:
        d=z[:,lo:hi]; b=u[lo:hi];row=rec[key]; a=eliminate(d)[0];ab=eliminate(np.vstack([d,b]))[0]
        assert (7000+a,7000+ab)==(row['rank'],row['augmented_rank']);i,j=row['quotient_pivot_rows'],row['quotient_pivot_cols'];n,det=eliminate(d[np.ix_(i,j)]);assert n==a and det==row['quotient_determinant'];v=cert[name]
        if name=='side_fit':
            v=v[rec['free_recipe_rows']]; assert not ((v@d.astype(np.int64)-b)%P).any();pair=None
        else:
            t=v[7000:];assert not (d.astype(np.int64)@t%P).any();pair=int(b.astype(np.int64)@t%P);assert pair==row['pairing']!=0
        report[name]={'rank':a,'augmented_rank':ab,'minor_determinant':det,'pairing':pair}
    return report

def direct(folder,wide,quotient):
    rec,cert=read(folder);H=array(folder/'independent10.npy',(7616,800));S=array(folder/'side.npy',(7616,312));y=array(folder/'independent_rhs.npy',(800,));b=array(wide/'rhs.npy',(7000,));u=cert['side_fit'];v=cert['fresh10_only'];t=cert['all_five'];lhs=np.zeros(8112,dtype=np.int64)
    assert 8112*(P-1)**2+P < 2**63-1
    hashes=[hashlib.sha256() for _ in range(2)];offset=0
    for n,path,digest in zip((6848,768),(wide/'matrix.npy',quotient/'new_rows.npy'),rec['wide_recipe_sha256']):
        X=array(path,(n,7000));h=hashlib.sha256()
        for k in range(0,n,37):
            x=X[k:k+37];h.update(x.tobytes());sl=slice(offset+k,offset+k+len(x));s=S[sl];f=H[sl];hashes[0].update(s.tobytes());hashes[1].update(f.tobytes());x=x.astype(np.int64);s=s.astype(np.int64);f=f.astype(np.int64)
            assert not ((x@v[:7000]+f@v[7000:])%P).any();assert not ((x@t[:7000]+s@t[7000:7312]+f@t[7312:])%P).any();lhs=(lhs+np.r_[u[sl]@x,u[sl]@s,u[sl]@f])%P
        assert h.hexdigest()==digest;offset+=n
    assert offset==7616 and [h.hexdigest() for h in hashes]==[rec['side_sha256'],rec['fresh_matrix_sha256']];assert (hv(b),hv(y))==(rec['wide_rhs_sha256'],rec['fresh_rhs_sha256'])
    assert np.array_equal(lhs[:7000],b) and not lhs[7000:7312].any();r=(y.astype(np.int64)-lhs[7312:])%P; ix=np.flatnonzero(r);i=int(ix[0]);first=[i,int(y[i]),int(lhs[7312+i]),int(r[i])];assert [len(ix),hv(r),first]==[rec['side_fit_holdout'][k] for k in ('nonzero','residual_sha256','first')]
    bv=int((b.astype(np.int64)@v[:7000]+y.astype(np.int64)@v[7000:])%P);bt=int((b.astype(np.int64)@t[:7000]+y.astype(np.int64)@t[7312:])%P);assert (bv,bt)==(rec['fresh10_only']['pairing'],rec['all_five']['pairing'])
    return {'modulus': P,'side_fit_verified_equations':7312,'side_fit_holdout_nonzero':len(ix),'first_holdout':[int(a) for a in first],'original_recipe_rows_annihilated_per_separator':[offset,offset],'fresh_only_pairing':bv,'all_five_pairing':bt,'side_sha256':hashes[0].hexdigest(),'fresh_sha256':hashes[1].hexdigest()}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--folder',type=Path,default=Path(__file__).resolve().parent/'four_source/five');p.add_argument('--wide',type=Path);p.add_argument('--quotient',type=Path);a=p.parse_args();assert bool(a.wide)==bool(a.quotient);print(json.dumps(direct(a.folder,a.wide,a.quotient) if a.wide else compact(a.folder),indent=2))
