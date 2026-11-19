"""Exact extension of the SU(3) j=34 profile-10 pivot minor to four coupled side profiles.

Every stage retains the original 7,000 profile-10 coordinates. Independent integer
products certify the selected frame only. The new 80-point profile-10 frame is
fixed in five_frames.py and is never used to choose the side-only coefficients.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import five_frames as frames
import profile10_missing as m
import four_source_refit as refit
f,fs,P=m.f,m.fs,m.P

def value(p):return json.loads(p.read_bytes())
def array(p):return np.load(p,mmap_mode='r',allow_pickle=False)
def digest(a):return fs.harray(a)

def archived(a,p,ww):
    rec=value(a.wide/'wide.json'); q=value(a.last/'result.json');b=array(a.wide/'rhs.npy');
    assert fs.fsha(a.wide/'wide.json')==q['old_record_sha256'] and fs.fsha(a.last/'result.json')==frames.QUOTIENT and fs.fsha(a.first/'build.json')==refit.FROZEN['build.json']
    old=array(a.wide/'matrix.npy');new=array(a.last/'new_rows.npy');stacks=[old,new];assert [digest(t) for t in (old,new,b)]==[rec['matrix_sha256'],q['new_rows_data_sha256'],rec['rhs_sha256']]
    rows=np.vstack(stacks);assert rows.shape==(7616,7000) and b.shape==(7000,)
    cols=frames.indices(p);ii=[k for t in frames.SIDE for k in cols[t]];z=cols[10]+ii
    e=a.archive/'edges.npy';k=a.archive/'kernel.npy';n=a.first/'new_rows.npy'; assert [fs.fsha(path) for path in (e,k,n)]==[f.joint.FILES['edges'],f.joint.FILES['kernel'],value(a.first/'build.json')['new_matrix_sha256']]
    edges=array(e);carry=array(k);other=array(n);assert (edges.shape,carry.shape,other.shape)==((2752,2244),(1519,2244),(4096,5044))
    assert not carry[:,z].any();assert not f.joint.rhs(a.root)[ii].any()
    # On arbitrary profile-10 points the carried directions are zero by source incidence.
    # v.inputs, called by frames.recipes, verified the complete carried source list.
    plus=frames.sides(a.runtime,p,ww[6848:],a.workers); checks=frames.direct_checks(a.runtime,p,ww[6848:],plus)
    side=np.vstack([edges[:,ii],other[:,[2800+j for j in ii]],plus]);assert side.shape==(7616,312)
    ci=q['quotient_pivot_cols'];F=q['free_coordinates'];oldcol=set(rec['pivot_columns']);assert F==[j for j in range(7000) if j not in oldcol] and sorted(ci)==list(range(365))
    I=rec['pivot_rows']+[6848+x for x in q['quotient_pivot_rows']];J=rec['pivot_columns']+[F[x] for x in ci];assert len(set(I))==len(set(J))==len(I)==len(J)==7000
    piv=dict(pivot_rows=I,pivot_columns=J,panel_orders=rec['panel_orders']+[365],panel_determinants=rec['panel_determinants']+[q['quotient_determinant']],minor_determinant=q['new_minor_determinant']);assert piv['minor_determinant']==102708
    meta=dict(prime=P,recipes=len(ww),carried_rows=1519,carried_array_sha256=f.joint.FILES['kernel'],carried_zero_entries_checked=len(z)*len(carry),carried_zero_counts={str(t):int(np.count_nonzero(carry[:,cols[t]])) for t in (*frames.SIDE,10)},side_targets=list(frames.SIDE),side_widths=[len(cols[t]) for t in frames.SIDE],side_sha256=digest(side),side_new_rows_sha256=digest(plus),side_explicit_checks=checks,wide_recipe_sha256=[digest(old),digest(new)],wide_rhs_sha256=digest(b),pivot_parent_receipts=[q['old_record_sha256'],frames.QUOTIENT],wide_minor_determinant=piv['minor_determinant'])
    return rows,b,side,piv,meta

def solve(X,b0,extra,side_rhs,chosen,factors,free_rows,Z,u,out,label):
    """Z has one row for each nonpivot recipe minus its projection on the pivot rows."""
    I,J,R,K,U,L,panels=factors; r=len(I); width=len(chosen);D=Z[:,chosen];y=u[chosen]; B=extra[:,chosen];b=side_rhs[chosen];K=K[:,chosen]
    i=f.joint.pivot_indices(f.matrix(D).transpose());j=f.joint.pivot_indices(f.matrix(D[i])) if i else []; n=len(i);assert n==len(j)==int(f.matrix(D).rank())
    rank1=int(f.matrix(np.vstack([D,y])).rank());assert rank1 in (n,n+1)
    M=f.matrix(D[np.ix_(i,j)]) if n else None; md=int(M.det()) if n else 1;assert md
    gamma=np.array([int(v[0]) for v in M.transpose().solve(f.matrix(y[j,None])).tolist()],dtype=np.int64) if n else np.zeros(0,dtype=np.int64)
    delta=(y.astype(np.int64)-gamma@D[i].astype(np.int64))%P;assert not delta[j].any() and bool(delta.any())==(n!=rank1)
    audit=dict(width=r+width,rank=r+n,augmented_rank=r+rank1,quotient_rank=n,quotient_augmented_rank=rank1,quotient_pivot_rows=i,quotient_pivot_cols=j,quotient_determinant=md,pivot_determinant=102708*md%P,quotient_rows=len(free_rows),quotient_columns=width)
    if rank1==n:
        cert=np.zeros(len(X),dtype=np.int64);cert[np.asarray(free_rows)[i]]=gamma;rhs=(b0.astype(np.int64)-refit.multiply(cert,X))%P;cert[I]=m.coefficients(rhs[J],U,L,panels)
        np.testing.assert_array_equal(refit.multiply(cert,X),b0);np.testing.assert_array_equal(refit.multiply(cert,B),b)
        audit.update(certificate='coefficients',coefficient_nonzeros=int(np.count_nonzero(cert)),integer_coordinates_checked=r+width,residual_nonzeros=0)
    else:
        h=int(np.flatnonzero(delta)[0]);t=np.zeros(width,dtype=np.int64);t[h]=1
        if n:t[j]=[-int(v[0])%P for v in M.solve(f.matrix(D[i,h,None])).tolist()]
        cert=np.zeros(r+width,dtype=np.int64);cert[r:]=t;cert[J]=(-K.astype(np.int64)@t)%P
        assert (r+width)*(P-1)**2<2**63
        for s in range(0,len(X),48):assert not ((X[s:s+48].astype(np.int64)@cert[:r]+B[s:s+48].astype(np.int64)@t)%P).any()
        pair=int((b0.astype(np.int64)@cert[:r]+b.astype(np.int64)@t)%P);assert pair==int(delta[h])!=0
        audit.update(certificate='separator',pairing=pair,independent_rows_annihilated=len(X),extra_column=h,augmented_minor_determinant=audit['pivot_determinant']*pair%P)
    cert=fs.packed(cert);audit['certificate_file_sha256']=fs.write_array(out/(label+'.npy'),cert)
    print(label,audit['rank'],audit['augmented_rank'],audit['certificate'],audit.get('pairing'),flush=True);return audit,cert

def main(a):
    sys.path.insert(0,str(a.runtime));p,ww=frames.recipes(a.root,a.first,a.last);a.out.mkdir(exist_ok=True,parents=True)
    X,b0,S,piv,meta=archived(a,p,ww);print('side frames and kernel checked',meta['side_sha256'],flush=True)
    assert not (a.out/'result.json').exists();fs.write_array(a.out/'side.npy',S)
    H,y=frames.fresh_frame(a.root,a.runtime,p,ww,a.workers);assert H.shape==(7616,800);fs.write_array(a.out/'independent10.npy',H);fs.write_array(a.out/'independent_rhs.npy',y)
    E=np.hstack((S,H));b=np.r_[np.zeros(S.shape[1],dtype='<u4'),y];A=np.hstack([X,E]);factors=m.kernel(A,piv);I,J,R,K,U,L,panels=factors;assert R==list(range(7000,7312+800))
    free=sorted(set(range(len(X)))-set(I));assert len(free)==616;Z=m.reduce(A[free],J,R,K);u=m.reduce(np.r_[b0,b][None,:],J,R,K)[0]
    for off in range(0,len(I),96):assert not m.reduce(A[I[off:off+96]],J,R,K).any()
    assert Z.shape==(616,1112);fs.write_array(a.out/'schur.npy',Z);fs.write_array(a.out/'schur_rhs.npy',u)
    meta.update(fresh_profile10_seed=frames.SEED,fresh_points=frames.POINTS,fresh_matrix_sha256=digest(H),fresh_rhs_sha256=digest(y),schur_sha256=digest(Z),schur_rhs_sha256=digest(u),free_recipe_rows=free)
    # Satisfy the four zero targets before opening the holdout in the fit.
    st,co=solve(X,b0,E,b,list(range(312)),factors,free,Z,u,a.out,'side_fit');meta['five_profiles']=st
    if st['certificate']=='coefficients':
        lhs=refit.multiply(co.astype(np.int64),H);err=(y.astype(np.int64)-lhs)%P;idx=np.flatnonzero(err);meta['side_fit_holdout']=dict(nonzero=len(idx),width=len(y),residual_sha256=digest(err),first=[int(idx[0]),int(y[idx[0]]),int(lhs[idx[0]]),int(err[idx[0]])] if len(idx) else None)
    # Compare without side equations to localize an obstruction on the widened frame.
    groups={'fresh10_only':list(range(312,1112)),'all_five':list(range(1112))}
    for k,sel in groups.items(): meta[k],_=solve(X,b0,E,b,sel,factors,free,Z,u,a.out,k)
    fs.putjson(a.out/'result.json',meta)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('root','first','last','wide','archive','runtime','out'):ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--workers',type=int,default=16);a=ap.parse_args();assert __debug__;main(a)
