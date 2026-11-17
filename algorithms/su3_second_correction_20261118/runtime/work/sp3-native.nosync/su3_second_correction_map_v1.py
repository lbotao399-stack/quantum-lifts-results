"""Q0 blocks in the existing second-target basis frames.

Preserve every original target letter matrix and add only missing source or
Q0-rule letters. A complete saved target basis then certifies the same frame.
"""
import sp_environment
import su_trace_adapter_v1 as algebra
from native_su_projection_v1 import project
from native_su_parallel_v1 import evaluate_profile
from g2_full_classical_screen_v1 import contents,sector
from fortuity.full.evalg import Frame
from fortuity.full.letters import Q_letter
from pathlib import Path
from math import ceil,comb
from functools import lru_cache
import argparse
import hashlib
import json
import pickle
import time
import numpy as np

HERE=Path(__file__).resolve().parent
PRIME=1000003


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


@lru_cache(maxsize=8)
def alphabets(mu,source_length,target_length):
    original=sorted({letter for content in contents(mu,target_length) for letter in content})
    full=set(original)|{letter for content in contents(mu,source_length) for letter in content}
    while True:
        following=full|{value for letter in full for _,left,right in Q_letter(letter) for value in (left,right)}
        if following==full:break
        full=following
    return original,sorted(full)


def target_frame(plan,index,prime=PRIME,points_override=0):
    target=plan['targets'][index];degree=target['fermions'];k=degree+2
    points=points_override or ceil((target['dimension']+32)/comb(k,degree))
    original,full=alphabets(tuple(plan['mu']),plan['source_length'],plan['target_length'])
    seed=2026111723+index*7919
    frame=Frame(3,points,k,[degree],original,seed=seed,p=prime)
    additional=[letter for letter in full if letter not in frame.L]
    extra_seed=2026111897+index*7919
    if additional:
        extra=Frame(3,points,k,[degree],additional,seed=extra_seed,p=prime)
        frame.L.update(extra.L)
    description=dict(N=3,points=points,k=k,degrees=[degree],alphabet=original,
        seed=seed,prime=prime,columns=frame.ncols,additional_letters=additional,additional_seed=extra_seed,
        preserved_target_basis_frame=True)
    return frame,description


def run(index,basis,output,start_row=0,stop_row=0,batch=8,seconds=900,workers=8):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=True)
    plan_path=HERE/'su3_second_correction_profiles_v1.json'
    plan=json.loads(plan_path.read_text());source=plan['sources'][index]
    source_path=basis/'independent_seeds.pkl';declared=json.loads((basis/'result.json').read_text())
    assert sha(source_path)==declared['seeds_sha256']
    seeds=pickle.loads(source_path.read_bytes());available=len(seeds)
    stop_row=stop_row or available
    assert 0<=start_row<stop_row<=available
    seeds=seeds[start_row:stop_row]
    recipes=pickle.dumps(seeds,protocol=4);recipe_path=output/'source_recipes.pkl'
    if recipe_path.exists():assert recipe_path.read_bytes()==recipes
    else:recipe_path.write_bytes(recipes)
    frames={j:target_frame(plan,j) for j in source['targets']}
    manifest=dict(source=source,prime=PRIME,mu=plan['mu'],length=plan['source_length'],parity=-1,
        source_recipes_sha256=sha(recipe_path),start_row=start_row,stop_row=stop_row,rows=len(seeds),
        complete_source_dimension=source['dimension'],targets={j:d for j,(_,d) in frames.items()},batch=batch,
        driver_sha256=sha(__file__),parallel_driver_sha256=sha(HERE/'native_su_parallel_v1.py'),
        plan_sha256=sha(plan_path))
    manifest=json.loads(json.dumps(manifest));marker=output/'manifest.json'
    if marker.exists():assert json.loads(marker.read_text())==manifest
    else:save(marker,manifest)
    records=[]
    for start in range(start_row,stop_row,batch):
        stop=min(start+batch,stop_row)
        missing=[j for j in source['targets'] if not (output/('t%02d_%05d_%05d.json'%(j,start,stop))).exists()]
        projection_seconds=0
        if missing:
            before=time.monotonic()
            expressions=[project(algebra.parity_project({word:1},-1),
                sector(tuple(plan['mu']),plan['source_length']),PRIME) for word in seeds[start-start_row:stop-start_row]]
            projection_seconds=time.monotonic()-before
        for j in source['targets']:
            name='t%02d_%05d_%05d'%(j,start,stop)
            marker,matrix_path=output/(name+'.json'),output/(name+'.npy')
            if marker.exists():
                record=json.loads(marker.read_text());assert sha(matrix_path)==record['sha256']
            else:
                frame=frames[j][0];before=time.monotonic()
                matrix=evaluate_profile(expressions,frame,source['profile'],plan['targets'][j]['profile'],workers)
                elapsed=time.monotonic()-before
                assert matrix.shape==(stop-start,frame.ncols)
                temp=matrix_path.with_suffix('.tmp.npy');np.save(temp,matrix.astype(np.uint32));temp.replace(matrix_path)
                record=dict(source=index,target=j,start=start,stop=stop,columns=frame.ncols,
                    projection_seconds=projection_seconds,evaluation_seconds=elapsed,
                    projected_terms=sum(map(len,expressions)),nonzero_entries=int(np.count_nonzero(matrix)),
                    workers=workers,sha256=sha(matrix_path))
                save(marker,record)
            records.append(record)
        save(output/'progress.json',dict(source=index,completed_blocks=len(records),
            required_blocks=ceil(len(seeds)/batch)*len(source['targets']),latest=records[-1],seconds=time.monotonic()-started))
        if time.monotonic()-started>=seconds:break
    required=ceil(len(seeds)/batch)*len(source['targets'])
    result=dict(status='complete_selected_evaluation_blocks' if len(records)==required else 'partial_evaluation_blocks',
        source=index,start_row=start_row,stop_row=stop_row,input_rows=len(seeds),complete_source_dimension=source['dimension'],
        completed_blocks=len(records),required_blocks=required,records=records,seconds=time.monotonic()-started,
        all_structural_targets=source['targets'],correction_solved=False,
        scope='Q0 blocks in preserved target-basis frames. No second-order solution or homogeneous quotient is inferred.')
    save(output/'result.json',result)
    print(json.dumps({key:value for key,value in result.items() if key!='records'}),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--index',type=int,required=True);parser.add_argument('--basis',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--start-row',type=int,default=0)
    parser.add_argument('--stop-row',type=int,default=0);parser.add_argument('--batch',type=int,default=8)
    parser.add_argument('--seconds',type=int,default=900);parser.add_argument('--workers',type=int,default=8)
    run(**vars(parser.parse_args()))
