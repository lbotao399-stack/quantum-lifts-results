"""Joint rank on identical source rows. Verify every discarded target column before concatenating.

Mode zero uses all eleven homogeneous equations and the source-18 action on
11 induced by the target-12 elimination. Other modes check independent subsets.
The reported rank always refers to the displayed row and column spaces only.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from flint import nmod_mat
import full_strand as fs
P=fs.P
SIZES=((1,501),(2,0),(4,1845),(10,1554),(6,225))


def column_factor(matrix, preferred=(), panel=80):
    """Select original independent columns, using the remaining rows after each block.

    Changes of pivot row are carried through an explicit row permutation. Each
    panel subtracts multiples of its normalized rows from every remaining row.
    With no remaining nonzero row, the original column rank is exactly known.
    """
    z=fs.packed(matrix).copy(); n,w=z.shape; assert z.size==0 or z.max()<P
    perm=np.arange(n);cols=[];used=set();preferred=list(dict.fromkeys(map(int,preferred)));assert all(0<=k<w for k in preferred)
    first=set(preferred);order=preferred+[i for i in range(w) if i not in first]
    r=0;ds=[];det=1
    for b in range(0,w,panel):
        if r==n:break
        batch=order[b:b+panel];block=z[r:,batch];
        if not block.any():continue
        R,h=nmod_mat(block.T.tolist(),P).rref();h=int(h);assert h
        select=[];u=0
        for j in range(h):
            while not R[j,u]:u+=1
            select.append(r+u);u+=1
        reduced,h2=nmod_mat(z[np.ix_(select,batch)].tolist(),P).rref();assert h2==h
        pc=[];u=0
        for j in range(h):
            while not reduced[j,u]:u+=1
            pc.append(batch[u]);u+=1
        # Bring selected rows to the front without copying the remaining large matrix.
        for j,src in enumerate(select):
            dst=r+j
            if src!=dst:
                z[[src,dst]]=z[[dst,src]];perm[[src,dst]]=perm[[dst,src]]
        stop=r+h;A=nmod_mat(z[r:stop,pc].tolist(),P);d=int(A.det());assert d
        U=fs.packed(fs.product(fs.packed(A.inv().tolist()),z[r:stop]));np.testing.assert_array_equal(U[:,pc],np.eye(h,dtype='<u4'));z[r:stop]=U
        for s in range(stop,n,256):
            v=z[s:s+256];v[:]=(v.astype(np.int64)-fs.product(v[:,pc].copy(),U))%P
        if stop<n:assert not z[stop:,batch].any()
        cols.extend(pc);used.update(pc);det=det*d%P;ds.append(d);r=stop
        if not np.any(z[r:]):break
    assert not z[r:].any() and len(cols)==r==len(used)
    rows=list(map(int,perm[:r]));return rows,cols,det,ds


def span_projection(arr, labels, out, t, hints=None):
    # Compress exactly, by proving every original column is in the span of
    # the retained ones. Independent original columns have the same left kernel.
    arr=np.asarray(arr,dtype='<u4');assert len(arr)==len(labels);n,w=arr.shape
    if hints is not None:
        rr=[labels.index(tuple(v)) for v in hints['rows']];cc=list(map(int,hints['columns']));d=hints['determinant']
        sub=arr[np.ix_(rr,cc)];assert fs.harray(sub)==hints['entries_sha256'] and len(cc)==len(rr)==hints['order']
    else:rr,cc,d,_=column_factor(arr,panel=64)
    I=arr[:,cc];piv=arr[rr,:];U,L,panels,actual,det=fs.factor(piv,cc,panel=80) if rr else (None,None,[],[],1)
    # A valid minor may have singular leading panels; changed coordinates are legal.
    if actual!=cc:
        exact=nmod_mat(arr[np.ix_(rr,cc)].tolist(),P);assert int(exact.det())==d;cc=list(actual);I=arr[:,cc];d=det
    assert det==d
    for b in range(0,n,96):
        v=arr[b:b+96]
        if not rr:assert not v.any();continue
        C,Y=fs.project(v,U,L,panels);assert not Y.any();np.testing.assert_array_equal(fs.product(C,piv),v)
    return I,[(t,j) for j in cc],dict(target=t,width=w,rank=len(cc),signed_minor_in_selected_order=d,selected_columns=cc,selected_rows=[labels[i] for i in rr],original_sha256=fs.harray(arr),all_coordinates_checked=n*w,compressed_sha256=fs.harray(I),source='strand quotient' if hints else 'side image')


def joint(blocks,keys,labels,out):
    from exact_blas_rows import matmul_mod
    Y=fs.packed(np.concatenate(blocks,axis=1));assert Y.shape==(len(labels),len(keys))
    rr,cc,d,pan=column_factor(Y);U,L,Ps,newcols,d1=fs.factor(Y[rr,:],cc,panel=80);assert len(newcols)==len(cc)
    if newcols!=cc:assert nmod_mat(Y[np.ix_(rr,cc)].tolist(),P).det()==d
    else:assert d==d1
    coeff=np.zeros((len(labels),len(rr)),dtype='<u4');sel=set(rr)
    for k in range(0,len(Y),100):
        v=Y[k:k+100];c,res=fs.project(v,U,L,Ps);assert not res.any();np.testing.assert_array_equal(matmul_mod(c,Y[rr,:],P),v);coeff[k:k+len(v)]=c
    np.testing.assert_array_equal(coeff[rr],np.eye(len(rr),dtype='<u4'))
    digest=fs.write_array(out/'shared_coefficients.npy',coeff);fre=[i for i in range(len(labels)) if i not in sel]
    return dict(rank=len(cc),source_rows=len(Y),quotient_parameter_count=len(fre),columns=len(keys),minor_determinant=d1,minor_sha256=fs.harray(Y[np.ix_(rr,newcols)]),pivot_source_rows=[labels[j] for j in rr],pivot_target_columns=[keys[j] for j in newcols],free_source_rows=[labels[j] for j in fre],coefficients_sha256=digest,coefficient_sign='free row minus its stored combination of pivot rows',compressed_entries_checked=int(Y.size),initial_elimination_pivots=d,panels=pan)


def validate(a):
    sys.path.insert(0,str(a.runtime.resolve()));from tempfile import TemporaryDirectory
    # Blocks with ranks 2 and 2 share a row direction; their joint rank is 3.
    with TemporaryDirectory() as tmp:
        path=Path(tmp);v=np.array([[1,2,0,3],[3,6,0,9],[0,0,1,2],[3,6,1,11],[0,0,0,0]],dtype='<u4');w=np.array([[1,0,3],[3,0,9],[0,0,0],[3,0,9],[0,1,7]],dtype='<u4');labels=[(0,i) for i in range(len(v))];v[[0,4]]=v[[4,0]];w[[0,4]]=w[[4,0]];b0,k0,r0=span_projection(v,labels,path,0);b1,k1,r1=span_projection(w,labels,path,1);res=joint([b0,b1],k0+k1,labels,path)
        assert r0['rank']==r1['rank']==2 and res['rank']==3 and res['quotient_parameter_count']==2;rr,cc,d,_=column_factor(v,[1,3,0],panel=1);assert nmod_mat(v[np.ix_(rr,cc)].tolist(),P).det()==d
        rng=np.random.default_rng(111254);checks=0
        for nr,nc,rk in ((9,13,4),(13,9,8),(7,11,7),(8,5,0),(5,8,1)):
            Z=fs.packed(fs.product(rng.integers(0,P,size=(nr,rk),dtype=np.int64),rng.integers(0,P,size=(rk,nc),dtype=np.int64)));exact=nmod_mat(Z.tolist(),P);assert exact.rank()==rk
            for panel in (1,2,4,7):
                ii,jj,val,_=column_factor(Z,list(reversed(range(nc))),panel);assert len(ii)==len(jj)==rk and (int(nmod_mat(Z[np.ix_(ii,jj)].tolist(),P).det()) if rk else 1)==val;checks+=1
    print(json.dumps(dict(individual=[2,2],joint=3,kernel=2,exact_reference_minors=checks)),flush=True)


def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));out=a.out;out.mkdir(exist_ok=True,parents=True);labels=[(s,i) for s,n in SIZES for i in range(n)];assert len(labels)==4125
    bs=[];keys=[];recs=[];ts=[7,8,9,12] if a.mode=='four' else [0,2,3,4,7,8,9,12,11]
    if a.mode=='zero':assert a.side11 is not None
    for t in ts:
        if t in fs.STRANDS:
            base=a.full/f'target_{t}';r=json.loads((base/'result.json').read_bytes());fp=base/'residual.npy';assert (r['target'],r['prime'])==(t,P) and r['residual_file_sha256']==fs.fsha(fp);A=np.load(fp,mmap_mode='r',allow_pickle=False);rows=labels if t==9 else labels[:3900];assert A.shape==(len(rows),r['full_width']);z,k,rec=span_projection(A,rows,out,t,r['quotient_minor'])
        else:
            base=(a.side11 if t==11 else a.reduced)/f'target_{t}';r=json.loads((base/'result.json').read_bytes());assert r['target']==t and r['prime']==P;arr=np.zeros((3900,r['width']),dtype='<u4');i=0
            if t==11:assert r['source18_carry_included']
            for b in r['rows']:
                s,n=SIZES[[s for s,n in SIZES].index(b['source'])];assert b['free']==n
                if not b['structural_zero']:
                    fp=base/b['path'];assert fs.fsha(fp)==b['sha256'];v=np.load(fp,mmap_mode='r',allow_pickle=False);assert v.shape==(n,r['width']);arr[i:i+n]=v
                i+=n
            assert i==3900;rows=labels[:3900];z,k,rec=span_projection(arr,rows,out,t)
        if t!=9:z=np.vstack([z,np.zeros((225,z.shape[1]),dtype='<u4')]);
        bs.append(z);keys.extend(k);recs.append(rec);fs.putjson(out/f'profile_{t}.json',rec);print(t,rec['rank'],flush=True)
    result=joint(bs,keys,labels,out);result.update(targets=ts,zero_by_previous_elimination=[1,5],prime=P,plan_sha256=fs.PLAN,mode=a.mode);fs.putjson(out/'result.json',result);print('joint',result['rank'],result['quotient_parameter_count'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['validate','four','zero']);ap.add_argument('--runtime',type=Path,default=fs.maps.RUNTIME);ap.add_argument('--full',type=Path);ap.add_argument('--reduced',type=Path);ap.add_argument('--side11',type=Path);ap.add_argument('--out',type=Path);a=ap.parse_args();assert __debug__;validate(a) if a.mode=='validate' else run(a)
