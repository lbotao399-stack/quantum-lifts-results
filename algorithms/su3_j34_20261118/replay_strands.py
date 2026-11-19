"""Replay a saved narrow strand minor in row batches without allocating a Python-list square matrix."""
import argparse,gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
from flint import nmod_mat
import projected_q0 as m
from kernel_strands import STRANDS

def replay(root,runtime,task,out,base):
    base=base or root/'strands'/('task_'+str(task)); record=json.loads((base/'result.json').read_bytes());blob=(base/'minor.json.gz').read_bytes();assert hashlib.sha256(blob).hexdigest()==record['sha256'];cert=json.loads(gzip.decompress(blob));pre,src,tgt=STRANDS[task%4];s,t=(pre,src) if task<4 else (src,tgt)
    data=json.loads((root/('classical_length9.json' if task<4 else 'classical_order1.json')).read_bytes());raw=(json.dumps({'order':1,'incoming':data},indent=2)+'\n').encode() if task<4 else (root/'classical_order1.json').read_bytes();plan=json.loads(raw)['incoming'];assert hashlib.sha256(raw).hexdigest()==record['plan_sha256'];assert cert['map']=={k:record[k] for k in cert['map']} and record['sources']==[s] and record['targets']==[t] and record['prime']==1000003
    m.setup(runtime.resolve(),plan,[t],[s],record['prime'],True);cols=cert['columns'];words=[(i,m.decode(w)) for i,w in cert['rows']];n=len(cols);assert n==len(set(cols))==len(set(words))==len(words)==record['rank_lower_bound'] and all(i==0 for i,w in words) and m.frames[0].ncols==record['frame_width']
    A=nmod_mat(n,n,record['prime']);hh=hashlib.sha256();
    for k in range(0,n,16):
        a=m.values(words[k:k+16])[:,cols];hh.update(a.astype('<u4').tobytes());
        for u,row in enumerate(a):
            for j,v in enumerate(row.tolist()):A[k+u,j]=v
    d=int(A.det());assert d==cert['determinant']==record['determinant'] and d!=0
    result=dict(task=task,sources=[s],targets=[t],prime=record['prime'],order=n,determinant=d,minor_sha256=record['sha256'],regenerated_entries=n*n,entries_le_u4_rowmajor_sha256=hh.hexdigest());out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=m.HERE);p.add_argument('--runtime',type=Path,default=m.RUNTIME);p.add_argument('--task',type=int,choices=range(8),required=True);p.add_argument('--base',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();replay(a.root,a.runtime,a.task,a.out,a.base)
