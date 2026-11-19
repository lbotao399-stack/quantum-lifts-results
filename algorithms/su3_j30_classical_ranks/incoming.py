"""Rank only the length-eight classical map in the spin-singlet R-007 second target.
Arguments are the existing p=1000033 production directory, published matrix-hash
receipt, and a new checkpoint directory. One column per source, one row per
physical evaluation coordinate. Never assemble Q1 or the inhomogeneous column.
"""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np

P = 1000033
PLAN = 'b8c669fa244e851f38a25b91186c3091e98a939935242805b0f71a9f8c00d609'
HASHES = 'ebac8b013044e79a3b7bf246555b2c68d3207f2798a580a8f0059b4c0845beb9'

def load(f): return json.loads(f.read_bytes())

def digest(f):
    h=hashlib.sha256()
    with f.open('rb') as g:
        for b in iter(lambda:g.read(4194304),b''): h.update(b)
    return h.hexdigest()

def save(f,x):
    tmp=f.with_suffix('.part');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(f)

def save_array(f,x):
    tmp=f.with_suffix('.part')
    with tmp.open('wb') as s: np.save(s,x,allow_pickle=False)
    tmp.replace(f)

def control():
    from bounded_flint_rref import reduce_target
    from flint import nmod_mat
    # An omitted source profile and a duplicate row across two targets. Compare
    # the embedded staged system with a direct dense rank over the same field.
    starts=[[0,1,4],[1,2,3,4],[],[2]]; matrices=[[[1,2,0],[0,1,3]],[[1,0,0,3],[2,1,1,0]],[],[[1]]]
    big=[];stages=[];n=5
    for cols,block in zip(starts,matrices):
        value=np.array(block,dtype=np.int32).reshape(len(block),len(cols));values=np.zeros((len(value),len(cols)+1),dtype=np.int32);values[:,:-1]=value
        full=np.zeros((len(value),n),dtype=np.int32);full[:,cols]=value;big.extend(full.tolist())
        st=reduce_target(values,np.array([*cols,n],dtype=np.int64),stages,n,P);assert st['consistent'];stages.append(st)
        assert int(nmod_mat(big,P).rank())==sum(x['rank'] for x in stages)
    return [x['rank'] for x in stages]

def run(root,receipt,out,threads):
    runtime=root/'runtime/work/sp3-native.nosync';sys.path.insert(0,str(runtime));import sp_environment
    from bounded_flint_rref import reduce_target
    from flint import ctx
    assert __debug__ and 1<=threads<=64;ctx.threads=threads
    test=control();plan_file=runtime/'su3_second_correction_profiles_v1.json'
    assert digest(plan_file)==PLAN and digest(receipt)==HASHES
    plan=load(plan_file);rec=load(receipt); assert rec['prime']==P and plan['mu']==[2,2,5,2,2] and plan['parity']==-1
    src=plan['sources'];tgt=plan['targets'];n=plan['source_dimension']; assert (n,plan['target_dimension'],len(src),len(tgt))==(49328,63983,57,37)
    offsets=np.cumsum([0]+[s['dimension'] for s in src]).tolist();assert offsets[-1]==n
    records={r['target']:r for r in rec['targets'] if r['equation']==2};assert set(records)==set(range(len(tgt)))
    stages=[];runs=[];stamps=[];out.mkdir(parents=True,exist_ok=True);packed=root/'output/compacted/second_Q0/maps'
    # Smallest first, exactly as the coupled solver does. Column indices always refer to the fixed complete source basis.
    for a in sorted(tgt,key=lambda t:(t['dimension'],t['index'])):
        j=a['index'];r=records[j];width=r['coordinates'];edges=[m for m in r['matrices'] if m[0]=='second_Q0'];expected=[s['index'] for s in src if j in s['targets']];assert [s[1] for s in edges]==expected
        marker=out/('stage%02d.json'%j);ar=out/('rows%02d.npy'%j);ac=out/('columns%02d.npy'%j)
        guard=dict(prime=P,target=j,input_digests=[e[-1] for e in edges],previous=stamps.copy(),plan_sha256=PLAN,driver_sha256=digest(Path(__file__)))
        if marker.exists():
            rs=load(marker);assert rs['guard']==guard and digest(ar)==rs['rows_sha256'] and digest(ac)==rs['columns_sha256']
        else:
            c=[v for i in expected for v in range(offsets[i],offsets[i+1])]+[n];val=np.zeros((width,len(c)),dtype=np.int32);pos=0
            for _,i,t,h,w,sha in edges:
                assert (t,h,w)==(j,src[i]['dimension'],width);file=packed/('profile%02d'%i)/('rows00000_%05d'%h)/('t%02d_00000_%05d.npy'%(j,h));assert digest(file)==sha
                x=np.load(file,allow_pickle=False,mmap_mode='r');assert x.dtype==np.dtype('<u4') and x.shape==(h,w) and x.max()<P
                val[:,pos:pos+h]=x.T;pos+=h
            assert pos==len(c)-1
            stage=reduce_target(val,np.asarray(c,dtype=np.int64),stages,n,P);assert stage['consistent'] and n not in stage['pivots'];del val
            save_array(ar,stage['rows']);save_array(ac,stage['columns']);rs=dict(guard=guard,pivots=stage['pivots'],rank=stage['rank'],rows_sha256=digest(ar),columns_sha256=digest(ac));save(marker,rs);del stage
        st=dict(rows=np.load(ar,allow_pickle=False,mmap_mode='r'),columns=np.load(ac,allow_pickle=False),pivots=rs['pivots'],rank=rs['rank']);assert st['rows'].shape==(st['rank'],len(st['columns'])) and not set(st['pivots']).intersection(p for s in stages for p in s['pivots'])
        stages.append(st);stamps.append(digest(marker));runs.append(dict(target=j,coordinates=width,matrices=len(edges),rank_increment=st['rank'],stage_sha256=stamps[-1]));save(out/'progress.json',dict(processed=len(runs),incoming_rank_lower=sum(q['rank_increment'] for q in runs),last=j))
    rank=sum(r['rank_increment'] for r in runs);assert len(runs)==37 and sum(r['matrices'] for r in runs)==221 and sum(r['coordinates'] for r in runs)==65406
    result=dict(prime=P,input_plan_sha256=PLAN,matrix_receipt_sha256=HASHES,driver_sha256=digest(Path(__file__)),control_increments=test,source_dimension=n,middle_dimension=63983,following_dimension=46876,stages=runs,incoming_rank_over_prime=rank,required_outgoing_plus_Cartan_rank_for_rational_exclusion=63983-rank,rational_incoming_rank_is_only_a_lower_bound=True);save(out/'result.json',result);print('Incoming classical rank',rank,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('production',type=Path);p.add_argument('receipt',type=Path);p.add_argument('output',type=Path);p.add_argument('--threads',type=int,default=1);a=p.parse_args();run(a.production.resolve(),a.receipt.resolve(),a.output.resolve(),a.threads)
