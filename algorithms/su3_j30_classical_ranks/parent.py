"""Extract a Cartan-vanishing length-nine cocycle and test it in the original incoming frames.

All coordinates are residues modulo 1000033. A rank jump certifies classical
non-exactness over this field only. The charge-sector rank and all input hashes
are independently pinned by the preceding public computations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

P = 1000033
N = 63983
M = 49328

def read(p): return json.loads(p.read_bytes())
def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def digest(v): return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def save(p,v):
    t=p.with_suffix('.tmp'); t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def array(p,v):
    t=p.with_suffix('.part')
    with t.open('wb') as f:np.save(f,v,allow_pickle=False)
    t.replace(p)
def loadmat(p,h,shape,dtype):
    assert sha(p)==h,p
    x=np.load(p,allow_pickle=False,mmap_mode='r');assert x.shape==shape and x.dtype==np.dtype(dtype),p
    return x

def dot(A,x):
    # Blockwise signed int64, reducing bounded dot products before accumulation.
    assert A.shape[1]==len(x) and x.dtype==np.int64 and np.all((x>=0)&(x<P))
    assert A.shape[1]*(P-1)**2<2**63
    out=np.empty(len(A),dtype=np.int64)
    for i in range(0,len(A),128):
        s=A[i:i+128]; assert not s.size or (s.min()>=0 and s.max()<P)
        out[i:i+128]=(s@x)%P
    return out

def kernel(stages,n,seed):
    piv={p for s in stages for p in s['pivots']};assert sum(len(s['pivots']) for s in stages)==len(piv) and all(0<=p<n for p in piv)
    free=np.array(sorted(set(range(n))-piv),dtype=np.int64);v=np.zeros(n+1,dtype=np.int64)
    v[free]=np.random.Generator(np.random.PCG64(seed)).integers(0,P,len(free),dtype=np.int64)
    for s in reversed(stages):
        ps=s['pivots']; assert not np.any(v[ps]) and s['columns'][-1]==n
        v[ps]=(-dot(s['rows'],v[s['columns']]))%P
    for s in stages: assert not dot(s['rows'],v[s['columns']]).any()
    return v[:-1],len(free)

def control():
    from bounded_flint_rref import reduce_target
    # Includes compacted non-contiguous columns, overlap between stages and rank zero.
    stages=[]
    for a,c in [([[1,2,0,0],[0,0,1,0]],[0,2,3,5]),([[1,1,1,0]],[0,1,4,5]),([[0,0]],[1,5])]:
        st=reduce_target(np.array(a,dtype=np.int32),np.array(c,dtype=np.int64),stages,5,P);stages.append(st)
    v,nfree=kernel(stages,5,11);assert nfree==2 and (v[0]+2*v[2])%P==v[3]==(v[0]+v[1]+v[4])%P==0
    from flint import nmod_mat
    assert int(nmod_mat([[1,0,0],[0,1,0],[2,1,1]],P).rank())==3
    for b in [0,1]:
        ss=reduce_target(np.array([[1,0,1],[0,1,1]],dtype=np.int32),np.arange(3),[],2,P)
        tt=reduce_target(np.array([[2,1,3+b]],dtype=np.int32),np.arange(3),[ss],2,P)
        assert tt['rank']==b and tt['consistent']==(b==0)
    return True

def run(base,pins,out,seed):
    jdir=base/'su3_j30_joint_rank_20261119/output';rows=base/'su3_j30_joint_map_20261119/output';idir=base/'su3_j30_incoming_rank_20261119/output';frames=base/'su3_j30_frames_fixed_20261119/output';root=base/'su3_prime1000033_v1_20261118';rt=root/'runtime/work/sp3-native.nosync';sys.path.insert(0,str(rt));import sp_environment
    from flint import ctx
    from bounded_flint_rref import reduce_target
    ctx.threads=64;assert __debug__ and control();pins=read(pins);out.mkdir(exist_ok=True,parents=True)
    joint=read(jdir/'result.json');incoming=read(idir/'result.json');assert joint['prime']==incoming['prime']==P and joint['stages']==pins['joint'] and incoming['stages']==pins['incoming']
    planpath=root/'runtime/work/sp3-native.nosync/su3_second_correction_profiles_v1.json';assert sha(planpath)==incoming['input_plan_sha256']==pins['original_plan_sha256'];plan=read(planpath)
    assert (plan['source_dimension'],plan['target_dimension'],len(plan['sources']),len(plan['targets']))==(M,N,57,37)
    receipt_file=idir.parent/'matrices.json';assert sha(receipt_file)==incoming['matrix_receipt_sha256']==pins['matrices_sha256'];edges={t['target']:t['matrices'] for t in read(receipt_file)['targets'] if t['equation']==2}
    offset=np.cumsum([0]+[s['dimension'] for s in plan['targets']]);soffset=np.cumsum([0]+[s['dimension'] for s in plan['sources']]);assert (offset[-1],soffset[-1])==(N,M)
    stages=[];markers={}
    for s in joint['stages']:
        k=s['stage'];f=jdir/(k+'.json');assert sha(f)==s['stage_sha256'];r=read(f);markers[k]=r
        cs=np.load(jdir/(k+'_columns.npy'),allow_pickle=False);assert sha(jdir/(k+'_columns.npy'))==r['columns_sha256'] and cs[-1]==N and np.all(cs[1:]>cs[:-1]);ps=r['pivots'];assert len(ps)==r['rank']==s['rank_increment'];a=loadmat(jdir/(k+'_rows.npy'),r['rows_sha256'],(len(ps),len(cs)),'<i4');stages.append(dict(pivots=ps,rows=a,columns=cs))
    v,nfree=kernel(stages,N,seed);assert nfree==32145;del stages;array(out/'parent.npy',v.astype('<u4'));print('parent extracted',len(v),nfree,flush=True)
    # Verify every raw outgoing and Cartan row chunk, and the ordered hashes seen by the reduction.
    raw={k:np.zeros(r['guard']['width'],dtype=np.int64) for k,r in markers.items()};inputs={k:[] for k in raw};counts={k:0 for k in raw}
    for i,t in enumerate(plan['targets']):
        assert i==t['index'];start=0
        for f in sorted((rows/f's{i:02d}').glob('*_*/result.json')):
            r=read(f);lo,hi=r['start'],r['stop'];assert lo==start<hi<=t['dimension'] and r['guard']['recipe']==pins['recipes'][i];start=hi
            for mk,m in r['matrices'].items():
                k=f'c{i:02d}' if mk=='c' else mk
                assert i in markers[k]['guard']['sources'] and m['shape']==[hi-lo,len(raw[k])];a=loadmat(f.parent/(mk+'.npy'),m['sha256'],tuple(m['shape']),'<u4')
                raw[k]=(raw[k]+dot(a.T,v[offset[i]+lo:offset[i]+hi]))%P;inputs[k].append([i,lo,hi,m['sha256']]);counts[k]+=1
        assert start==t['dimension']
    for k,r in markers.items():assert digest(inputs[k])==r['guard']['blocks'] and not np.any(raw[k]),k
    print('all outgoing and Cartan residuals zero',flush=True)
    # Independently evaluate the chosen source in each of the 37 original incoming evaluation frames.
    vals={}; frame_rows=0
    for i,t in enumerate(plan['targets']):
        f=frames/f'e2_t{i:02d}';mf=f/'manifest.json';rf=f/'result.json';assert [sha(mf),sha(rf)]==pins['frames'][i];meta=read(mf);rec=read(rf)
        assert meta['recipe_sha256']==pins['recipes'][i] and meta['target']==t and meta['prime']==rec['prime']==P and rec['manifest_sha256']==sha(mf) and rec['seed_rank']==rec['dimension']==t['dimension'] and rec['injective']
        val=np.zeros(rec['width'],dtype=np.int64);start=0
        for ck in rec['chunks']:
            lo,hi=ck['start'],ck['stop'];assert lo==start<hi<=t['dimension'];start=hi
            a=loadmat(f/f'rows{lo:05d}_{hi:05d}.npy',ck['sha256'],(hi-lo,len(val)),'<u4');val=(val+dot(a.T,v[offset[i]+lo:offset[i]+hi]))%P
        assert start==t['dimension'];vals[i]=val;frame_rows+=start
    evaluated=np.concatenate([vals[i] for i in range(37)]).astype('<u4');assert frame_rows==N and len(evaluated)==65406 and np.any(evaluated);array(out/'evaluated.npy',evaluated)
    print('all 37 incoming frames evaluated',flush=True)
    report=dict(prime=P,seed=seed,free=nfree,joint_rank=N-nfree,incoming_rank=incoming['incoming_rank_over_prime'],parent_dimension=N,parent_nonzero=int(np.count_nonzero(v)),recipe_vector_sha256=sha(out/'parent.npy'),evaluation_sha256=sha(out/'evaluated.npy'),evaluation_nonzero=int(np.count_nonzero(evaluated)),evaluations=len(evaluated),outgoing_and_Cartan=[dict(key=k,coordinates=len(raw[k]),chunks=counts[k],nonzero=int(np.count_nonzero(raw[k]))) for k in raw],incoming_checks=[])
    assert report['incoming_rank']==32144;save(out/'closure.json',report);stages=[];jump=False
    for rec in incoming['stages']:
        i=rec['target'];mk=idir/f'stage{i:02d}.json';assert sha(mk)==rec['stage_sha256'];old=read(mk);es=[e for e in edges[i] if e[0]=='second_Q0'];ids=[e[1] for e in es];assert ids==[s['index'] for s in plan['sources'] if i in s['targets']] and len(es)==rec['matrices'] and old['guard']['input_digests']==[e[-1] for e in es]
        cols=[j for k in ids for j in range(soffset[k],soffset[k+1])]+[M];w=len(vals[i]);assert rec['coordinates']==w
        x=np.zeros((w,len(cols)),dtype=np.int32);pos=0
        for _,k,ix,n,c,h in es:
            assert (ix,c,n)==(i,w,plan['sources'][k]['dimension']);f=root/'output/compacted/second_Q0/maps'/f'profile{k:02d}'/f'rows00000_{n:05d}'/f't{i:02d}_00000_{n:05d}.npy';a=loadmat(f,h,(n,w),'<u4');x[:,pos:pos+n]=a.T;pos+=n
        assert pos==len(cols)-1;x[:,-1]=vals[i];s=reduce_target(x,np.asarray(cols,dtype=np.int64),stages,M,P);del x
        assert [p for p in s['pivots'] if p!=M]==old['pivots'] and old['rank']==rec['rank_increment'];jump=not s['consistent'];assert s['rank']==old['rank']+int(jump)
        if jump:assert s['pivots'][-1]==M and int(np.count_nonzero(s['rows'][-1]))==1 and int(s['rows'][-1,-1])==1
        report['incoming_checks'].append(dict(target=i,old_increment=old['rank'],augmented_increment=s['rank'],inconsistent=jump));save(out/'prefix.json',report);print('target',i,'augmented jump',int(jump),flush=True)
        if jump:break
        stages.append(s)
    # Inconsistency in any prefix proves inconsistency in the full system; appending one column adds at most one rank.
    assert jump or len(stages)==37
    report.update(nonboundary_mod_prime=jump,augmented_rank=report['incoming_rank']+int(jump),quantum_lifts=0);save(out/'result.json',report)
    print('full augmented rank',report['augmented_rank'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('base',type=Path);p.add_argument('pins',type=Path);p.add_argument('out',type=Path);p.add_argument('--seed',type=int,default=2026111937);a=p.parse_args();run(a.base,a.pins,a.out,a.seed)
