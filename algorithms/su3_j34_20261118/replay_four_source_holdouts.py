"""Regenerate both held-out profile-10 equations from the sparse frozen coefficients.

Use the original, unsliced frame and accumulate projected source recipes directly.
No training matrix, elimination output other than the coefficients, or Slurm file is needed.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import four_source_joint as joint
import four_source_refit as refit


def replay(root, runtime, which):
    sys.path.insert(0,str(runtime))
    import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from native_su_values_v1 import evaluate
    v, fs, prime = joint.v, joint.fs, joint.P
    path = root/'four_source'
    plan, old, _ = v.inputs(root, root/'allcoord')
    assert fs.fsha(path/'new_words.json.gz') == refit.FROZEN['new_words.json.gz']
    words = old + joint.read_words(path/'new_words.json.gz'); assert len(set(words)) == len(words) == 6848
    terms = list(v.rhs_terms(root, plan).items()); records=[]
    for name, directory, rank_file, expected_file, seed in (
        ('first',path,'fit.json','validation.json',joint.VALIDATION_SEED),
        ('second',path/'refit','rank.json','fresh.json',refit.FRESH_SEED)):
        if which not in (name,'both'):continue
        rank=json.loads((directory/rank_file).read_bytes()); expect=json.loads((directory/expected_file).read_bytes())
        cpath=directory/'coefficients.npy'; csha=fs.fsha(cpath)
        assert csha==rank['certificate_file_sha256']==expect.get('frozen_coefficients_sha256',expect.get('frozen_solution_sha256'))
        assert rank.get('precommitted_fresh_seed',rank.get('precommitted_validation_seed'))==expect['seed']==seed
        c=np.load(cpath,allow_pickle=False).astype(np.int64); assert len(c)==8367 and c.min()>=0 and c.max()<prime
        cc=np.r_[c[:2752],c[4271:]]; fr,charge=v.frame(plan,expect['points'],seed); target=v.maps.pattern(plan['targets'][10]); rhs=np.zeros(fr.ncols,dtype=np.int64);lhs=rhs.copy()
        for batch in v.batches(terms,176): rhs=(rhs+evaluate([dict(batch)],fr)[0])%prime
        for s,indices,seq in joint.tasks_for(words,11):
            keep=[(idx,w) for idx,w in zip(indices,seq) if cc[idx]]
            if not keep: continue
            ex=[project(alg.parity_project({w:1},1),charge,prime) for _,w in keep]; coeffs=cc[[i for i,_ in keep]]
            lhs=(lhs + coeffs@derivative(ex,fr,v.maps.pattern(plan['sources'][s]),target))%prime
        residual=(rhs-lhs)%prime
        for key,val in [('rhs',rhs),('lhs',lhs),('residual',residual)]: np.testing.assert_array_equal(val,expect[key])
        assert expect['width']==len(lhs)==120 and expect['nonzero']==int(np.count_nonzero(residual))==120 and fs.harray(residual)==expect['residual_sha256']
        records.append(dict(fit=name,prime=prime,seed=seed,coordinates=len(lhs),nonzero=int(np.count_nonzero(residual)),terms=len(terms),nonzero_profile10_coefficients=int(np.count_nonzero(cc)),coefficients_sha256=csha,rhs_sha256=fs.harray(rhs),lhs_sha256=fs.harray(lhs),residual_sha256=fs.harray(residual)))
    return records

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--which',choices=('both','first','second'),default='both');ap.add_argument('--output',type=Path);args=ap.parse_args(); assert __debug__
    root=Path(__file__).resolve().parent; rows=replay(root,joint.v.maps.RUNTIME,args.which)
    if args.output:joint.fs.putjson(args.output,rows)
    print(json.dumps(rows,indent=2))
