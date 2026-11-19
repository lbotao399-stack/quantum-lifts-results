"""Pair new, rationally projected source directions with the archived modular separator."""
import argparse, hashlib, json, random, sys
from pathlib import Path
import numpy as np
import profile10_span as span
from profile10_validation import EDGE

SHA='ae565ce3f91869bac4dfbb026dda9b185e23591490af4995970bbecd557814d1'
SEED=761590321

def compute(root,rt):
    v=span.v; cert=root/'allcoord'; f=cert/'profile10_span.json'; assert span.fsha(f)==SHA
    rec=json.loads(f.read_bytes());p=rec['prime']; assert rec['annihilated_recipe_rows']==2752 and rec['pairing']==779716
    span.setup(rt,root,cert,rec['points'],rec['seed']); assert rec['shape']==[len(span.rows),span.fr.ncols]
    from g2_full_classical_screen_v1 import contents
    from su3_stream_rank_v3 import random_word
    from sp_full_profiles import profile
    from fortuity.full.evalg import nferm
    rng=random.Random(SEED);used=set(span.rows);c=np.asarray(rec['annihilator'],dtype=np.int64)
    assert len(c)==span.fr.ncols and len(c)*(p-1)**2<2**63
    cases=[];plan=span.plan
    for s in EDGE:
        pat=v.maps.pattern(plan['sources'][s]); choices=[x for x in contents(tuple(plan['mu']),plan['source_length']) if profile(x)==pat]
        for attempt in range(1,1001):
            w=random_word(rng.choice(choices),rng,8)
            if not w or (s,w) in used:continue
            e=span.project(span.alg.parity_project({w:1},1),span.charge,p)
            if not e:continue
            vals=np.concatenate([span.derivative([e],f,pat,span.tpat)[0] for f in span.frames]);pair=int(vals@c%p)
            if not pair:continue
            explicit={u:int(k)%p for u,k in span.alg.differential(e).items() if k%p and profile(tuple(x for t in u for x in t))==span.tpat}
            # Independent Taylor differential and Python exterior-algebra evaluation on every coordinate.
            words=sorted(explicit); z=np.zeros(span.fr.ncols,dtype=np.int64)
            for group in v.batches(words,24):
                co=np.asarray([explicit[u] for u in group],dtype=np.int64)
                for j,f in enumerate(span.frames):
                    i=16*j*span.fr.cpp;z[i:i+f.ncols]=(z[i:i+f.ncols]+co@f.eval_words(group,dtype=np.int64))%p
            np.testing.assert_array_equal(z,vals);assert pair==int(z@c%p)
            cases.append(dict(source_profile=s,source_dimension=plan['sources'][s]['dimension'],attempt=attempt,recipe=w,projected_terms=len(e),target_terms=len(explicit),evaluations_sha256=span.harray(span.packed(z)),pairing=pair));print(s,pair,attempt,flush=True);break
        else:raise ValueError('No witness found for source '+str(s))
    return dict(prime=p,certificate_sha256=SHA,frame_seed=rec['seed'],points=rec['points'],coordinates=len(c),sampling_seed=SEED,verified_by='profile derivation and explicit differential with Python exterior evaluation on all coordinates',archived_rhs_pairing=rec['pairing'],cases=cases)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=span.v.maps.HERE);ap.add_argument('--runtime',type=Path,default=span.v.maps.RUNTIME);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();assert __debug__;span.putjson(a.output,compute(a.root,a.runtime.resolve()))
