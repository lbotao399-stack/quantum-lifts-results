"""Bounded profile bases for a supplied full nonabelian SU(3) correction plan."""
import sp_environment
import su_trace_adapter_v1 as algebra
from su3_stream_rank_v3 import random_word
from su_modular_rowspace_v1 import RowSpace
from native_su_projection_v1 import project
from native_su_parallel_v1 import evaluate
from g2_full_classical_screen_v1 import contents, sector
from sp_full_profiles import profile
from fortuity.full.evalg import Frame, nferm
from fortuity.full.letters import letters_up_to
from pathlib import Path
from math import ceil, comb
import argparse
import hashlib
import json
import pickle
import random
import time
import numpy as np

HERE=Path(__file__).resolve().parent
PRIME=1000003


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(index,output,plan_path,role,seconds=300,batch=16,resume_from=None,workers=8):
    started=time.monotonic();output=Path(output);output.mkdir(parents=True,exist_ok=True)
    plan=json.loads(plan_path.read_text())
    MU=tuple(plan['mu']);length=plan['source_length' if role=='sources' else 'target_length']
    task=plan[role][index];assert task['index']==index
    selected_profile=tuple((tuple(key),n) for key,n in task['profile'])
    choices=[content for content in contents(MU,length) if profile(content)==selected_profile]
    assert choices
    degree,dimension=task['fermions'],task['dimension']
    k=degree+2
    alphabet=sorted({L for c in contents(MU,length) for L in c})
    frame=Frame(3,ceil((dimension+32)/comb(k,degree)),k,[degree],alphabet,
                seed=2026111723+index*7919,p=PRIME)
    manifest=dict(task=task,prime=PRIME,mu=MU,length=length,parity=-1,role=role,columns=frame.ncols,
                  frame=dict(N=3,points=frame.P,k=k,degrees=[degree],alphabet=alphabet,seed=2026111723+index*7919),
                  native_point_workers=workers,pending_term_limit=32768,parallel_driver_sha256=sha(HERE/'native_su_parallel_v1.py'),max_trace_length=6,rowspace_sha256=sha(HERE/'su_modular_rowspace_v1.py'),blas_sha256=sha(HERE/'exact_blas_rows.py'),driver_sha256=sha(__file__),plan_sha256=sha(plan_path))
    manifest=json.loads(json.dumps(manifest))
    marker=output/'manifest.json'
    if marker.exists():assert json.loads(marker.read_text())==manifest
    else:marker.write_text(json.dumps(manifest,indent=2)+'\n')
    space=RowSpace(frame.ncols,PRIME);seeds=[];seen=set();attempts=0;prior_seconds=0
    rng=random.Random(2026111741+index*1009)
    state_path=output/'checkpoint.pkl'
    loaded_path=state_path
    if not loaded_path.exists() and resume_from is not None and (resume_from/'checkpoint.pkl').exists():
        previous=json.loads((resume_from/'manifest.json').read_text())
        assert all(previous[k]==manifest[k] for k in ('task','prime','mu','length','parity','columns','frame','max_trace_length'))
        loaded_path=resume_from/'checkpoint.pkl'
        (output/'resumed_from.json').write_text(json.dumps(dict(path=str(loaded_path),sha256=sha(loaded_path)),indent=2)+'\n')
    if loaded_path.exists():
        old=pickle.loads(loaded_path.read_bytes());space=RowSpace.restore(old['space'])
        seeds,seen,attempts,prior_seconds=old['seeds'],old['seen'],old['attempts'],old['seconds']
        rng.setstate(old['rng'])
    pending=[];pending_seeds=[];pending_terms=0;peak_terms=0;last=started
    timings=dict(projection=0.0,evaluation=0.0,row_reduction=0.0)

    def flush():
        nonlocal pending_terms
        if not pending:return
        before=time.monotonic();values=evaluate(pending,frame,workers);timings['evaluation']+=time.monotonic()-before
        before=time.monotonic();chosen=space.extend(values);timings['row_reduction']+=time.monotonic()-before
        seeds.extend(pending_seeds[i] for i in chosen)
        assert len(seeds)==len(space.pivots)<=dimension
        pending.clear();pending_seeds.clear();pending_terms=0

    def checkpoint(final=False):
        state=dict(space=space.dump(),seeds=seeds,seen=seen,attempts=attempts,rng=rng.getstate(),
                   seconds=prior_seconds+time.monotonic()-started)
        temp=state_path.with_suffix('.tmp');temp.write_bytes(pickle.dumps(state,protocol=4));temp.replace(state_path)
        recipes=output/'independent_seeds.pkl';recipes_temp=recipes.with_suffix('.tmp');recipes_temp.write_bytes(pickle.dumps(seeds,protocol=4));recipes_temp.replace(recipes)
        result=dict(status='complete_rational_profile_basis_recipes' if len(seeds)==dimension else 'partial_profile_basis',
            profile_index=index,fermions=degree,dimension=dimension,rank_lower_bound=len(seeds),columns=frame.ncols,
            attempts=attempts,peak_projected_terms=peak_terms,timings=timings,seconds=state['seconds'],
            stopped_by='dimension_saturated' if len(seeds)==dimension else ('time_cap' if final else 'running'),
            seeds_sha256=sha(recipes),checkpoint_sha256=sha(state_path),driver_sha256=sha(__file__),
            scope='Independent modular evaluations certify independent rational projected trace recipes. Completeness is asserted only when the exact character dimension is reached. No Q0 image or correction is computed here.')
        temp=output/'result.json.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(output/'result.json')
        print(json.dumps(result),flush=True)

    while time.monotonic()-started<seconds and len(seeds)<dimension:
        attempts+=1
        word=random_word(rng.choice(choices),rng,6)
        if word is None or word in seen:continue
        seen.add(word)
        before=time.monotonic();expression=project(algebra.parity_project({word:1},-1),sector(MU,length),PRIME)
        timings['projection']+=time.monotonic()-before
        if not expression:continue
        assert all(profile(tuple(L for trace in w for L in trace))==selected_profile for w in expression)
        peak_terms=max(peak_terms,len(expression))
        pending.append(expression);pending_seeds.append(word);pending_terms+=len(expression)
        if len(pending)>=batch or pending_terms>=32768:flush()
        if time.monotonic()-last>30:
            flush();checkpoint();last=time.monotonic()
    flush();checkpoint(True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--index',type=int,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--seconds',type=int,default=300)
    ap.add_argument('--plan-path',type=Path,required=True);ap.add_argument('--role',choices=('sources','targets'),required=True)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--batch',type=int,default=16);ap.add_argument('--resume-from',type=Path);args=ap.parse_args()
    run(**vars(args))
