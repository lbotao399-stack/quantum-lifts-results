#!/usr/bin/env python3
"""Exact formal-trace versus finite-colour preparation for Sp(3) omega."""
from pathlib import Path
import sys,os,time,json,gzip,pickle,hashlib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
LOCAL_UP=HERE/'inputs/algorithms/dual_singlets_20261118'
DEFAULT_UP=ROOT/'quantum_lifts_frontier_work/root_audit/so11_pin/algorithms/dual_singlets_20261118'
UP=Path(os.environ.get('N4_DUAL_SOURCE',str(LOCAL_UP if (LOCAL_UP/'mixed_incoming.py').exists() else DEFAULT_UP)))
LOCAL_RT=HERE/'inputs/algorithms/su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
DEFAULT_RT=ROOT/'quantum_lifts_frontier_source/algorithms/su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
os.environ.setdefault('N4_TRACE_RUNTIME',str(LOCAL_RT if LOCAL_RT.exists() else DEFAULT_RT))
sys.path.insert(0,str(UP))
import numpy as np
from scipy.sparse import csr_matrix,save_npz
from flint import nmod_mat
import mixed_incoming as m
from trace_descent import highest_weights
from mixed33_exact import exact,apply

def mod_minor(rows):
    words=sorted({w for row in rows for w in row});ids={w:i for i,w in enumerate(words)}
    a=nmod_mat(len(rows),len(words),m.P)
    for i,row in enumerate(rows):
        for w,n in row.items():a[i,ids[w]]=n
    rr,rank=a.rref();piv=[];j=0
    for i in range(rank):
        while not rr[i,j]:j+=1
        piv.append(j);j+=1
    indep=m.pivotcols(a.transpose())
    sq=nmod_mat([[a[i,j] for j in piv] for i in indep],m.P)
    det=int(sq.det());assert det
    return {'rank':rank,'shape':[len(rows),len(words)],'basis_rows':indep,'columns':piv,'determinant':det,'words':words}

def main():
    start=time.monotonic();cache=HERE/'stable_rows.pkl'
    if cache.exists():
        parents,centre,incoming,outgoing=pickle.loads(cache.read_bytes())
    else:
        w=(3,3,5,0,0)
        parents=highest_weights(w,8)[3];centre=highest_weights(w,9)[4]
        print('formal dimensions',len(parents),len(centre),flush=True)
        incoming=[m.q(r) for r in parents];outgoing=[m.q(r) for r in centre]
        assert all(not m.q(r) for r in incoming+outgoing)
        cache.write_bytes(pickle.dumps((parents,centre,incoming,outgoing),protocol=4))
    print('rows ready',len(parents),len(centre),'seconds',time.monotonic()-start,flush=True)
    result={'weight':[3,3,5,0,0],'length':9,'prime':m.P}
    for name,rows in [('incoming',incoming),('outgoing',outgoing)]:
        info=mod_minor(rows)
        result[name]={k:v for k,v in info.items() if k!='words'}
        print(name,result[name]['shape'],info['rank'],info['determinant'],'seconds',time.monotonic()-start,flush=True)
    result['formal_centre_dimension']=len(centre)
    result['formal_cohomology_upper_bound']=len(centre)-result['incoming']['rank']-result['outgoing']['rank']
    result['nilpotence_checked_rows']=len(incoming)+len(outgoing)
    result['seconds']=time.monotonic()-start
    words=sorted({w for row in outgoing for w in row});ids={w:i for i,w in enumerate(words)}
    rr=[];cc=[];vv=[]
    for i,row in enumerate(outgoing):
        for w,n in row.items():rr.append(ids[w]);cc.append(i);vv.append(n)
    matrix=csr_matrix((np.asarray(vv,dtype=np.int64),(rr,cc)),shape=(len(words),len(centre)))
    save_npz(HERE/'stable_q_columns.npz',matrix)
    result['stable_Q0_sparse_shape']=list(matrix.shape)
    result['stable_Q0_sparse_nonzeros']=int(matrix.nnz)
    (HERE/'stable_rank_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('STABLE H upper',result['formal_cohomology_upper_bound'],flush=True)

if __name__=='__main__':main()
