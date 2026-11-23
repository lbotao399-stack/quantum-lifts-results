#!/usr/bin/env python3
"""Complete finite-colour relation space and its classical image for omega."""
from pathlib import Path
import json,time,hashlib
import numpy as np
from scipy.sparse import load_npz
from flint import nmod_mat

HERE=Path(__file__).resolve().parent
P=1000003

def pivots(a):
    r,rank=a.rref();out=[];j=0
    for i in range(rank):
        while not r[i,j]:j+=1
        out.append(j);j+=1
    return out

def main():
    start=time.monotonic()
    data=np.load(HERE/'finite_evaluation.npz')
    e=data['evaluation'];basis=data['basis_rows'];free=np.asarray(sorted(set(range(e.shape[0]))-set(map(int,basis))),dtype=np.int32)
    square=nmod_mat(e[basis].T.tolist(),P)
    coefficients=np.asarray(square.solve(nmod_mat(e[free].T.tolist(),P)).tolist(),dtype=np.int64).T
    relations=np.zeros((len(free),len(e)),dtype=np.int64)
    relations[np.arange(len(free)),free]=1
    relations[:,basis]=-coefficients%P
    print('relation solve',len(free),'seconds',time.monotonic()-start,flush=True)
    assert not any((nmod_mat(relations.tolist(),P)*nmod_mat(e.tolist(),P)).entries())
    assert np.array_equal(relations[:,free],np.eye(len(free),dtype=np.int64))
    q=load_npz(HERE/'stable_q_columns.npz')
    images=np.asarray(q@relations.T).T%P
    a=nmod_mat(images.tolist(),P)
    cols=pivots(a);rows=pivots(a.transpose());assert len(cols)==len(rows)
    minor=nmod_mat(images[np.ix_(rows,cols)].tolist(),P)
    determinant=int(minor.det());assert determinant
    finite=np.load(HERE/'finite_outgoing.npz')
    omega_lift=np.zeros(len(e),dtype=np.int64)
    omega_lift[basis]=finite['cycle']
    omega_defect=np.asarray(q@omega_lift)%P
    np.savez_compressed(HERE/'relation_image.npz',relations=relations.astype('uint32'),images=images.astype('uint32'),free_rows=free,basis_rows=basis,minor_rows=np.asarray(rows,dtype=np.int32),minor_columns=np.asarray(cols,dtype=np.int32),omega_lift=omega_lift.astype('uint32'),omega_defect=omega_defect.astype('uint32'))
    stable=json.loads((HERE/'stable_rank_result.json').read_text())
    native=json.loads((HERE/'native_frames_result.json').read_text())
    result={'passed':True,'prime':P,'stable_dimension':int(e.shape[0]),'finite_dimension':int(e.shape[1]),'complete_relation_dimension':len(free),'relation_kernel_residual_nonzero':0,'relation_image_rank':len(cols),'relation_image_minor_determinant':determinant,'stable_Q0_rank':stable['outgoing']['rank'],'finite_Q0_rank':native['outgoing_minor_rank'],'fortuitous_dimension':stable['outgoing']['rank']-len(cols)-native['outgoing_minor_rank'],'seconds':time.monotonic()-start,'scope':'All181 finite-field identities and all6020 formal output coordinates. Rational dimension follows separately from exact stable acyclicity and the archived rational finite class; finite-field relation coefficients are not asserted to be rational lifts.'}
    assert result['fortuitous_dimension']==1
    (HERE/'relation_image_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
