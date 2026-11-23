#!/usr/bin/env python3
"""Regenerate complete finite-colour evaluation and both adjacent maps."""
from probe_stable import HERE,m,np,nmod_mat
from pathlib import Path
import pickle,json,gzip,time,hashlib
from dual_charge_bands import coordinates,dense_frame

def main():
    start=time.monotonic()
    parents,centre,incoming,outgoing=pickle.loads((HERE/'stable_rows.pkl').read_bytes())
    page=json.loads(gzip.decompress((HERE/'inputs/algorithms/dual_singlets_20261118/sp3_omega/target.json.gz').read_bytes()))
    w=tuple(page['weight']);lib=m.library(HERE,6,m.P)
    letters=sorted({l for content in m.contents(w,9) for l in content})
    probes=page['centre_frames']['4']['samples']
    t=time.monotonic()
    data=coordinates(centre+incoming,4,letters,probes,6,'Sp',lib,1)
    e,b=data[:len(centre)],data[len(centre):]
    np.savez_compressed(HERE/'finite_evaluation.npz',evaluation=e.astype('uint32'),incoming=b.astype('uint32'),basis_rows=np.asarray(page['centre_frames']['4']['minor_rows'],dtype=np.int32))
    print('evaluated source+incoming',data.shape,'seconds',time.monotonic()-t,flush=True)
    basis=page['centre_frames']['4']['minor_rows']
    det=int(nmod_mat(e[basis].tolist(),m.P).det())
    assert det==page['centre_frames']['4']['determinant_mod_prime']
    ib=page['incoming'];bminor=b[np.ix_(ib['rows'],ib['columns'])]
    assert int(nmod_mat(bminor.tolist(),m.P).det())==ib['determinant']
    assert hashlib.sha256(bminor.astype('<u4').tobytes()).hexdigest()==ib['minor_sha256']
    functional=np.asarray(page['representative']['target_functional'],dtype=np.int64)
    assert not np.any(b@functional %m.P)
    cycle=np.asarray(page['representative']['cycle_basis_coefficients'],dtype=np.int64)
    assert int((cycle@e[basis]%m.P)@functional %m.P)==1
    print('finite minor and all boundary pairings passed',flush=True)
    nextletters=sorted({l for content in m.contents(w,10) for l in content})
    following=page['following_frames']['5']['samples']
    a=coordinates([outgoing[i] for i in basis],5,nextletters,following,6,'Sp',lib,1)
    oa=page['outgoing'];aminor=a[np.ix_(oa['rows'],oa['columns'])]
    assert int(nmod_mat(aminor.tolist(),m.P).det())==oa['determinant']
    assert hashlib.sha256(aminor.astype('<u4').tobytes()).hexdigest()==oa['minor_sha256']
    assert not np.any(cycle@a %m.P)
    np.savez_compressed(HERE/'finite_outgoing.npz',outgoing=a.astype('uint32'),cycle=cycle.astype('uint32'),functional=functional.astype('uint32'))
    # Separate dense Python trace evaluator on a reproducible coordinate subset.
    words=sorted({v for r in centre+incoming for v in r})
    selected=[words[i] for i in np.linspace(0,len(words)-1,12,dtype=int)]
    points,seed,keep=probes[0]
    f=dense_frame(6,'Sp',points,4,letters,seed)
    native=m.values(selected,f,lib)
    python=f.eval_words(selected,dtype=np.int64)
    np.testing.assert_array_equal(native,python)
    result={'passed':True,'prime':m.P,'formal_source_dimension':len(centre),'finite_source_dimension':len(basis),'evaluation_shape':list(e.shape),'all_boundary_rows':len(incoming),'boundary_pairings_nonzero':0,'omega_pairing':1,'outgoing_shape':list(a.shape),'outgoing_minor_rank':oa['rank'],'outgoing_minor_determinant':oa['determinant'],'incoming_minor_rank':ib['rank'],'incoming_minor_determinant':ib['determinant'],'finite_evaluation_determinant':det,'independent_dense_python_coordinates':int(native.size),'seconds':time.monotonic()-start}
    (HERE/'native_frames_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
