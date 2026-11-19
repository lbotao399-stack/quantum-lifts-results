"""Exact membership on all 12 evaluation points for the 4,271 fixed correction directions.

Uses only NumPy and python-flint. The two ranks are recomputed separately;
a positive membership saves coefficients, a negative one saves an annihilator.
The raw matrices are the two outputs of expanded_joint.py, in that order.
"""
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
from flint import nmod_mat

P=1000003
FILES={'edges':'e3c96352bede6ddd3f8c7501b0c7bbc1ba249a1b12ea68cc3fc02f5fe0f92468',
       'kernel':'c63270ea44213f133bf57539ef471e65308ea4e21ca2ea5f29bad268e9dfd373'}
PLAN='1b7ac3dec77fcb50473d01176a8f309773e8411f0fd4e969f58299a026df6e8c'
RHS='cb2cafc82eb5eb462d410520fef0c87ee026a221203dfbaa3c5f7f2f44c57330'

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def digest(v):return hashlib.sha256(np.ascontiguousarray(v,dtype='<u4').tobytes()).hexdigest()
def rhs(root):
    plan=root/'classical_order1.json';assert sha(plan)==PLAN
    z=json.loads(plan.read_bytes())['incoming']; b=[]
    for i,t in enumerate(z['targets']):
        assert i==t['index'];w=12*math.comb(t['fermions']+2,t['fermions'])
        if i in (6,10,13,14,15,16):
            path=root/('joint/edges' if i in (6,10) else 'first_profile16/coupled')/f'fresh_{i}.npy'
            record=json.loads(path.with_suffix('.json').read_bytes());assert sha(path)==record['vectors_sha256'] and record['prime']==P and record['seed']==88304773+229*i
            v=np.load(path,allow_pickle=False);v=v if i in (6,10) else v[2];assert v.shape==(w,) and digest(v)==record['residual_sha256']
        else:
            assert str(i) not in z['first_rhs_trace_terms'] and all(i not in z['sources'][s]['targets'] for s in (19,20,21,22))
            v=np.zeros(w,dtype='<u4')
        b.extend(v)
    b=np.array(b,dtype='<u4');assert digest(b)==RHS;return b

def matrix(a):return nmod_mat(a.tolist(),P)
def pivot_indices(m):
    z,r=m.rref();p=[];j=0
    for i in range(r):
        while not z[i,j]:j+=1
        assert int(z[i,j])==1;p.append(j);j+=1
    return p

def run(folder,root):
    b=rhs(root);blocks=[]
    for k,n in [('edges',2752),('kernel',1519)]:
        path=folder/(k+'.npy'); assert sha(path)==FILES[k];v=np.load(path,allow_pickle=False);assert v.shape==(n,len(b)) and v.dtype==np.dtype('<u4') and v.max()<P;blocks.append(v)
    A=np.vstack(blocks); del blocks; rows,width=A.shape; rank=int(matrix(A).rank()); aug=int(matrix(np.vstack([A,b])).rank())
    I=pivot_indices(matrix(A).transpose());assert len(I)==rank
    sel=A[I];J=pivot_indices(matrix(sel));assert len(J)==rank
    T=matrix(sel[:,J]);det=int(T.det());assert det
    coeff=np.array([int(x[0]) for x in T.transpose().solve(matrix(b[J,None])).tolist()],dtype=np.int64)
    # Largest inner dimension is at most 4271. Integer products cannot overflow int64.
    assert rows*(P-1)**2 < np.iinfo(np.int64).max
    residual=(b.astype(np.int64)-coeff@sel.astype(np.int64))%P
    assert not residual[J].any() and aug==rank+int(residual.any());bad=np.flatnonzero(residual)
    data=dict(prime=P,source_matrix_file_sha256=FILES,rhs_data_sha256=RHS,shape=[rows,width],rank=rank,augmented_rank=aug,pivot_rows=I,pivot_columns=J,pivot_minor_determinant=det,pivot_minor_data_sha256=digest(sel[:,J]),coefficients_on_pivot_rows=coeff.tolist(),coefficient_nonzeros=int(np.count_nonzero(coeff)),residual_nonzeros=len(bad),residual_sha256=digest(residual),membership=(aug==rank))
    if len(bad):
        j=int(bad[0]);dual=np.zeros(width,dtype=np.int64);dual[j]=1;dual[J]=[-int(x[0]) for x in T.solve(matrix(sel[:,j,None])).tolist()];dual%=P
        assert all(not ((v.astype(np.int64)@dual)%P).any() for v in np.array_split(A,24));v=int(b@dual%P);assert v==int(residual[j])!=0
        data.update(extra_column=j,annihilator=dual.tolist(),annihilator_checks=rows,pairing=v,augmented_minor_determinant=det*v%P)
    return data

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('folder',type=Path);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert __debug__;r=run(a.folder,a.root);a.output.write_text(json.dumps(r,separators=(',',':'))+'\n');print(r['rank'],r['augmented_rank'],r['coefficient_nonzeros'],r['residual_nonzeros'],flush=True)
