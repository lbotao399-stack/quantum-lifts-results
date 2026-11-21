"""Test selected incoming Q0 equations on exactly the same rows and columns as the descendant.

A separator proves nonexactness at the recorded prime. A common preimage proves only
consistency on the selected target profiles. Requires numpy and python-flint.
"""
import argparse
import hashlib
from pathlib import Path
import numpy as np
from flint import nmod_mat
from boundary import P, SELECTED, RECEIPT, get, plan, save, sha


def read_array(path, dim, pin):
    assert sha(path)==pin, path
    x=np.load(path,allow_pickle=False)
    assert x.dtype==np.dtype('<u4') and x.shape==dim and (x<P).all(),path
    return x


def mat(x):return nmod_mat(x.tolist(),P)


def pivots(x):
    if not x.shape[0] or not x.shape[1]:return []
    m,r=mat(x).rref(inplace=True);col=0;out=[]
    for i in range(r):
        while int(m[i,col])==0:col+=1
        assert int(m[i,col])==1;out.append(col);col+=1
    return out


def mul(x,y):
    # Every intermediate (including every partial sum) is an integer less than 2^53.
    assert x.shape[1]*(P-1)**2<2**53
    return ((x.astype(np.float64)@y.astype(np.float64)).astype(np.int64)%P).astype('<u4')


def solve(a,b):return np.array([[int(e) for e in row] for row in mat(a).solve(mat(b)).tolist()],dtype='<u4')


def analyse(a,b,source_order):
    assert b.shape==(1,a.shape[1]);initial=[];lo=0
    for s in source_order:
        hi=lo+s['dimension'];block=a[lo:hi];cols=np.flatnonzero(block.any(axis=0));initial.extend(lo+i for i in pivots(block[:,cols].T));lo=hi
    assert lo==len(a)
    rows=[initial[i] for i in pivots(a[initial].T)];cols=pivots(a[rows]);r=len(rows);assert r==len(cols)>0
    basis=a[rows];minor=basis[:,cols];det=int(mat(minor).det());assert det
    inv=np.array([[int(e) for e in row] for row in mat(minor).inv().tolist()],dtype='<u4')
    # Exact rank upper bound: every original row reconstructs, on ALL columns.
    for lo in range(0,len(a),64):
        block=a[lo:lo+64];assert np.array_equal(mul(mul(block[:,cols],inv),basis),block)
    c=mul(b[:,cols],inv);res=(b.astype(np.int64)-mul(c,basis))%P;out=dict(coordinates=a.shape[1],source_rows=len(a),rank=r,minor_determinant=det,minor_rows=rows,minor_columns=cols)
    nz=np.flatnonzero(res[0]);out['augmented_rank']=r+int(bool(len(nz)))
    if len(nz):
        j=int(nz[0]);f=np.zeros((a.shape[1],1),dtype='<u4');f[j,0]=1;f[cols,0]=(-mul(inv,basis[:,j:j+1]).astype(np.int64)[:,0])%P
        assert not mul(a,f).any() and int(mul(b,f)[0,0])==int(res[0,j])!=0
        out.update(verdict='nonexact_mod_prime',separator=f[:,0].tolist(),pairing=int(res[0,j]))
    else:
        pre=np.zeros((1,len(a)),dtype='<u4');pre[0,rows]=c[0];assert np.array_equal(mul(pre,a),b)
        out.update(verdict='consistent_on_selected_profiles',solution_on_minor_rows=c[0].tolist())
    return out


def assemble(op,ids,blocks,rhs,receipt,remote=None,extra_plan=None):
    a,_=plan(op);from boundary import PLANS
    assert sha(receipt)==RECEIPT and get(receipt)['descendants'][op]['plan_sha256']==PLANS[op]
    scope=[(r['index'],r['dimension']) for r in a['sources'] if set(ids)&set(r['targets'])];rhsrows=[];frames={};targetdata=[]
    def metadata(role,i):
        k=f'{op}/{role}{i:02d}';folder=blocks/k;x=get(folder/'done.json');s=a[role][i];m=get(folder/'manifest.json');r=get(folder/'result.json');assert (x['op'],x['index'],x['role'],x['dimension'],x['prime'],x['plan_sha256'])==(op,i,role,s['dimension'],P,PLANS[op]);assert m['task']==s and m['role']==role and m['prime']==P and m['plan_sha256']==PLANS[op] and r['rank_lower_bound']==r['dimension']==s['dimension'] and r['seeds_sha256']==x['seeds_sha256']==sha(folder/'independent_seeds.pkl');
        if remote is not None:assert x==remote[k],k
        return x,folder
    for j in ids:
        t=a['targets'][j];x,_=metadata('targets',j);f=rhs/f'{op}_t{j:02d}.npy';rr=get(f.with_suffix('.json'));raw=get(receipt)['descendants'][op]['profile_blocks'][j];assert (rr['op'],rr['prime'],rr['target'],rr['plan_sha256'],rr['input_sha256'],rr['terms'])==(op,P,j,PLANS[op],raw['sha256'],raw['terms']);assert rr['frame']==x['frame'];b=read_array(f,(1,rr['columns']),rr['values_sha256']);rhsrows.append(b);frames[j]=(x['frame'],rr['columns']);targetdata.append(dict(index=j,dimension=t['dimension'],columns=rr['columns'],frame=x['frame'],values_sha256=rr['values_sha256']))
    allrows=[];sources=[]
    for i,dim in scope:
        x,folder=metadata('sources',i);bl=[];tgs={};pred=a['sources'][i]
        for j in ids:
            tag,w=frames[j];has=i in a['targets'][j]['sources'];assert has==(j in pred['targets']);
            if not has:bl.append(np.zeros((dim,w),dtype='<u4'));continue
            m=x['targets'][str(j)];assert (tag,w)==(m['frame'],m['columns']);bl.append(read_array(folder/f't{j:02d}.npy',(dim,w),m['sha256']));tgs[j]=m['sha256']
        allrows.append(np.hstack(bl));sources.append(dict(index=i,dimension=dim,seeds_sha256=x['seeds_sha256'],blocks=tgs))
    arows=np.vstack(allrows);bb=np.hstack(rhsrows);answer=analyse(arows,bb,sources)
    answer.update(target_profiles=targetdata,sources=sources,prime=P,plan_sha256=PLANS[op],target_dimension=sum(r['dimension'] for r in targetdata),matrix_row_major_sha256=hashlib.sha256(arows.astype('<u4',copy=False).tobytes()).hexdigest(),rhs_row_major_sha256=hashlib.sha256(bb.astype('<u4',copy=False).tobytes()).hexdigest());return answer

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blocks',type=Path);p.add_argument('rhs',type=Path);p.add_argument('receipt',type=Path);p.add_argument('output',type=Path);p.add_argument('--a8-manifest',type=Path);args=p.parse_args();assert __debug__;remote=get(args.a8_manifest)['rows'] if args.a8_manifest else None;result={}
    for op,ids in SELECTED.items():
        r=assemble(op,ids,args.blocks,args.rhs,args.receipt,remote);result[op]=r;print(op,*[r[k] for k in ('source_rows','coordinates','rank','augmented_rank','verdict')],flush=True)
    save(args.output,dict(results=result,lowering_receipt_sha256=RECEIPT))
