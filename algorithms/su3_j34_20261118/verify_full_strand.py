"""Test a forced change of pivot columns against a direct finite-field inverse.
Optionally replay three archived SU(3) source words on an actual full target frame.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from flint import nmod_mat
import full_strand as fs


def verify(runtime):
    sys.path.insert(0, str(Path(runtime).resolve()))
    p = fs.P
    x = np.array([[0,2,4,1,0,6,2,1],[0,3,6,0,1,2,5,3],
                  [1,3,0,4,3,9,0,1],[0,1,0,2,2,5,7,4],
                  [8,3,13,19,27,31,4,7]], dtype='<u4')
    preferred = [0,2,1,3,4]
    assert nmod_mat(x[:2,preferred[:2]].tolist(),p).det()==0
    # The first five rows paired against themselves prescribe the identity
    # coefficients; the last three include two independent transverse directions.
    rhs=np.vstack([x,[2,5,8,1,0,3,4,9],[7,11,0,4,9,0,12,3],[1,1,2,0,7,9,3,15]]).astype('<u4')
    U,L,panels,cols,det=fs.factor(x,preferred,panel=2)
    assert cols[:2]!=preferred[:2] and L[2:,:2].any()
    X=nmod_mat(x.tolist(),p);B=nmod_mat(rhs.tolist(),p);M=nmod_mat(x[:,cols].tolist(),p)
    C=nmod_mat(rhs[:,cols].tolist(),p)*M.inv();R=B-C*X
    c,r=fs.project(rhs,U,L,panels);np.testing.assert_array_equal(c,C.tolist());np.testing.assert_array_equal(r,R.tolist());assert det==int(M.det());np.testing.assert_array_equal(c[:5],np.eye(5));assert not r[:,cols].any()
    rank=fs.quotient_rank(r,[j for j in range(x.shape[1]) if j not in cols],list(range(len(rhs))));assert rank['order']==int(R.rank());assert nmod_mat(r[np.ix_(rank['rows'],rank['columns'])].tolist(),p).det()==rank['determinant']
    return dict(prime=p,source_rows=5,target_coordinates=8,preferred=preferred,selected=cols,determinant=det,checked_coefficient_entries=int(c.size),checked_residual_entries=int(r.size),quotient_rank=rank['order'])


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--runtime',type=Path,default=fs.maps.RUNTIME);ap.add_argument('--actual',type=Path,help='Working directory for three full-frame source-15 trace rows');a=ap.parse_args();result=verify(a.runtime)
    if a.actual:
        root=Path(__file__).resolve().parent;f=root/'classical_order1.json';assert fs.fsha(f)==fs.PLAN;plan=json.loads(f.read_bytes())['incoming'];t=9;task,s,sh,_=fs.STRANDS[t];c,words=fs.readwords(root/f'strands/task_{task}/minor.json.gz',sh);a.actual.mkdir(parents=True,exist_ok=True)
        sel=[0,len(words)//2,len(words)-1];x,receipt=fs.build(a.actual/'sample.npy',[words[i] for i in sel],plan,t,s,a.runtime,1)
        from math import comb
        degree=plan['targets'][t]['fermions']; per_point=comb(degree+2,degree); points=(plan['targets'][t]['dimension']+8+per_point-1)//per_point
        assert x.shape==(3,points*per_point)==(3,3510) and receipt['direct_word_checks']==[[0,0,per_point],[1,points//2,per_point],[2,points-1,per_point]]
        result.update(actual_target=t,original_recipe_indices=sel,actual=receipt)
    print(json.dumps(result,indent=2))
