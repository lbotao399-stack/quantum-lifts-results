"""Reduce outgoing blocks first, then Cartan blocks, keeping the same 63,983 source indices."""
import argparse
import json
import math
import sys
from pathlib import Path
import numpy as np
import joint as j
import produce_joint


def width(t,cartan=False):
    d=t['fermions'];stride=math.comb(d+2,d);return stride*math.ceil((t['cartan_dimension' if cartan else 'dimension']+32)/stride)

def stages_for(plan):
    result=[]
    for t in sorted(plan['targets'],key=lambda t:(width(t),t['index'])):
        i=t['index']; result.append((f'q{i:02d}',[s['index'] for s in plan['sources'] if i in s['targets']],width(t)))
    for s in sorted((s for s in plan['sources'] if s['cartan_dimension']),key=lambda s:(width(s,True),s['index'])):
        result.append(('c%02d'%s['index'],[s['index']],width(s,True)))
    return result

def run(root,rows,out,threads):
    import incoming as inc
    from bounded_flint_rref import reduce_target
    from flint import ctx
    assert 1<=threads<=64;ctx.threads=threads;assert inc.control()==[2,1,0,1];j.inputs(root);maps=j.runtime(root);plan=j.load_plan();tasks=produce_joint.intervals(plan);sources=plan['sources'];offsets=np.cumsum([0]+[s['dimension'] for s in sources]);n=int(offsets[-1]);assert n==63983
    manifests={};guards={};intervals={i:[] for i in range(len(sources))};out.mkdir(exist_ok=True,parents=True)
    for i in range(len(sources)):
        fs=j.frames(maps,plan,i);guards[i]=j.make_guard(plan,i,fs);del fs
    for i,lo,hi in tasks:
        f=rows/f's{i:02d}'/f'{lo:05d}_{hi:05d}';r=json.loads((f/'result.json').read_bytes());assert r['guard']==guards[i] and (r['start'],r['stop'])==(lo,hi) and set(r['matrices'])==set(guards[i]['frames']);manifests[i,lo,hi]=r['matrices'];intervals[i].append((lo,hi))
    for s in sources:
        seq=sorted(intervals[s['index']]);assert seq[0][0]==0 and seq[-1][1]==s['dimension'] and all(a[1]==b[0] for a,b in zip(seq,seq[1:]))
    stages=[];receipts=[];stamps=[];outgoing=0
    for key,idx,w in stages_for(plan):
        mk=key if key[0]=='q' else 'c';files=[];ncol=0;cols=[]
        for i in idx:
            ncol+=sources[i]['dimension'];cols.extend(range(offsets[i],offsets[i+1]));
            for lo,hi in sorted(intervals[i]):
                f=rows/f's{i:02d}'/f'{lo:05d}_{hi:05d}'/(mk+'.npy');m=manifests[i,lo,hi][mk];assert m['shape']==[hi-lo,w] and j.sha(f)==m['sha256'];files.append((i,lo,hi,f,m['sha256']))
        cols.append(n);guard=dict(driver=j.sha(__file__),inputs=j.INPUTS,plan=j.PLAN,worker=j.sha(j.__file__),key=key,blocks=j.digest([[i,lo,hi,h] for i,lo,hi,_,h in files]),previous=stamps.copy(),width=w,sources=idx)
        ref=out/(key+'.json');rf=out/(key+'_rows.npy');cf=out/(key+'_columns.npy')
        if ref.exists():
            rec=json.loads(ref.read_bytes());assert rec['guard']==guard and j.sha(rf)==rec['rows_sha256'] and j.sha(cf)==rec['columns_sha256']
        else:
            v=np.zeros((w,ncol+1),dtype=np.int32);starts={s:k for s,k in zip(idx,np.cumsum([0]+[sources[i]['dimension'] for i in idx[:-1]]))}
            for i,lo,hi,f,h in files:
                x=np.load(f,allow_pickle=False,mmap_mode='r');assert x.dtype==np.dtype('<u4') and x.shape==(hi-lo,w) and x.max()<j.P;v[:,starts[i]+lo:starts[i]+hi]=x.T
            st=reduce_target(v,np.asarray(cols,dtype=np.int64),stages,n,j.P);assert st['consistent'] and n not in st['pivots'];del v
            inc.save_array(rf,st['rows']);inc.save_array(cf,st['columns']);rec=dict(guard=guard,rank=st['rank'],pivots=st['pivots'],rows_sha256=j.sha(rf),columns_sha256=j.sha(cf));j.save(ref,rec);del st
        st=dict(rows=np.load(rf,allow_pickle=False,mmap_mode='r'),columns=np.load(cf,allow_pickle=False),pivots=rec['pivots'],rank=rec['rank']);assert st['rows'].shape==(st['rank'],len(st['columns'])) and not set(st['pivots']).intersection(p for stage in stages for p in stage['pivots']);stages.append(st);stamps.append(j.sha(ref));receipts.append(dict(stage=key,coordinates=w,rank_increment=st['rank'],stage_sha256=stamps[-1]));
        if key[0]=='q':outgoing+=st['rank']
        joint=sum(s['rank'] for s in stages);assert joint<=31839
        j.save(out/'progress.json',dict(finished=len(stages),rank=joint,latest=receipts[-1])); print(key,st['rank'],joint,flush=True)
    assert len(stages)==51 and sum(r['rank_increment'] for r in receipts[:19])==outgoing
    result=dict(prime=j.P,plan=j.PLAN,physical_inputs=j.INPUTS,source_dimension=n,incoming_modular_lower_bound=32144,outgoing_measured_rank=outgoing,joint_measured_rank=joint,increment_from_Cartan=joint-outgoing,rational_cohomology_upper_bound=n-32144-outgoing,rational_Cartan_invisible_cohomology_upper_bound=n-32144-joint,intervals=len(tasks),Q0_blocks=sum(len(s['targets']) for s in sources),Cartan_blocks=sum(bool(s['cartan_dimension']) for s in sources),stages=receipts,quantum_lifts_certified_by_this_test=0);j.save(out/'result.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('row_dir',type=Path);p.add_argument('output',type=Path);p.add_argument('--threads',type=int,default=64);a=p.parse_args();sys.path.insert(0,str(a.root.resolve()/j.RT));import sp_environment;print(json.dumps(run(a.root,a.row_dir,a.output,a.threads)))
