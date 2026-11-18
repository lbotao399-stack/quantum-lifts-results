"""Nonabelian Q1 blocks for the complete first-correction coordinate space."""
from fractions import Fraction
from pathlib import Path
import argparse,json,pickle,time
import numpy as np
import su3_second_correction_map_v1 as maps
import su_quantum_integer_v1 as quantum
from native_su_parallel_v1 import evaluate
from sp_full_profiles import profile

HERE=maps.HERE
PRIME=maps.PRIME


def quantum_image(word,mu,prime=PRIME):
    projected=maps.project(maps.algebra.parity_project({word:1},-1),maps.sector(mu,9),prime)
    image=quantum.act(projected,1,3)
    modular={w:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime for w,c in image.items()}
    return projected,{w:c for w,c in modular.items() if c}


def run(index,basis,output,start_row=0,stop_row=0,batch=2,seconds=900,workers=8,points_override=0):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=True)
    plan_path=HERE/'su3_second_quantum_profiles_v1.json';plan=json.loads(plan_path.read_text())
    second=json.loads((HERE/'su3_second_correction_profiles_v1.json').read_text())
    source=plan['sources'][index]
    record=json.loads((basis/'result.json').read_text());seed_path=basis/'independent_seeds.pkl'
    assert record['dimension']==record['rank_lower_bound']==source['dimension']
    assert maps.sha(seed_path)==record['seeds_sha256']
    seeds=pickle.loads(seed_path.read_bytes());stop_row=stop_row or len(seeds)
    assert 0<=start_row<stop_row<=len(seeds)
    source_recipes=output/'source_recipes.pkl';data=pickle.dumps(seeds[start_row:stop_row],protocol=4)
    if source_recipes.exists():assert source_recipes.read_bytes()==data
    else:source_recipes.write_bytes(data)
    frames={j:maps.target_frame(second,j,points_override=points_override) for j in source['targets']}
    by_profile={tuple((tuple(key),n) for key,n in plan['targets'][j]['profile']):j for j in source['targets']}
    manifest=dict(source=source,prime=PRIME,source_mu=plan['source_mu'],target_mu=plan['target_mu'],
        start_row=start_row,stop_row=stop_row,batch=batch,points_override=points_override,
        source_recipes_sha256=maps.sha(source_recipes),source_basis_sha256=record['seeds_sha256'],
        targets={j:description for j,(_,description) in frames.items()},plan_sha256=maps.sha(plan_path),
        driver_sha256=maps.sha(__file__),quantum_driver_sha256=maps.sha(HERE/'su_quantum_integer_v1.py'),
        frame_driver_sha256=maps.sha(HERE/'su3_second_correction_map_v1.py'))
    manifest=json.loads(json.dumps(manifest));marker=output/'manifest.json'
    if marker.exists():assert json.loads(marker.read_text())==manifest
    else:maps.save(marker,manifest)
    records=[]
    for start in range(start_row,stop_row,batch):
        stop=min(start+batch,stop_row);expressions={j:[] for j in frames};before=time.monotonic()
        missing=[j for j in frames if not (output/('t%02d_%05d_%05d.json'%(j,start,stop))).exists()]
        quantum_seconds=0;input_terms=output_terms=0
        if missing:
            for word in seeds[start:stop]:
                projected,image=quantum_image(word,tuple(plan['source_mu']))
                input_terms+=len(projected);output_terms+=len(image)
                parts={j:{} for j in frames}
                for w,c in image.items():
                    key=profile(tuple(letter for trace in w for letter in trace))
                    assert key in by_profile,(index,key)
                    parts[by_profile[key]][w]=c
                for j in frames:expressions[j].append(parts[j])
            quantum_seconds=time.monotonic()-before
        for j,(frame,_) in frames.items():
            name='t%02d_%05d_%05d'%(j,start,stop);marker=output/(name+'.json');path=output/(name+'.npy')
            if marker.exists():
                result=json.loads(marker.read_text());assert maps.sha(path)==result['sha256']
            else:
                before=time.monotonic();values=evaluate(expressions[j],frame,workers)
                temp=path.with_suffix('.tmp.npy');np.save(temp,values.astype(np.uint32));temp.replace(path)
                result=dict(source=index,target=j,start=start,stop=stop,columns=frame.ncols,
                    quantum_seconds=quantum_seconds,evaluation_seconds=time.monotonic()-before,
                    projected_terms=input_terms,quantum_terms=output_terms,
                    target_terms=sum(map(len,expressions[j])),nonzero_entries=int(np.count_nonzero(values)),sha256=maps.sha(path))
                maps.save(marker,result)
            records.append(result)
        maps.save(output/'progress.json',dict(completed_rows=stop-start_row,required_rows=stop_row-start_row,
            completed_blocks=len(records),seconds=time.monotonic()-started))
        if time.monotonic()-started>=seconds:break
    required=((stop_row-start_row+batch-1)//batch)*len(frames)
    result=dict(status='complete_selected_Q1_evaluation_blocks' if len(records)==required else 'partial_Q1_evaluation_blocks',
        source=index,start_row=start_row,stop_row=stop_row,completed_blocks=len(records),required_blocks=required,
        records=records,seconds=time.monotonic()-started,points_override=points_override,
        scope='Q1 matrix blocks on all first-correction coordinates; no homogeneous parameter has been discarded. The coupled Q0 and Q1 equations remain to be solved.')
    maps.save(output/'result.json',result)
    print(json.dumps({key:value for key,value in result.items() if key!='records'}),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--index',type=int,required=True);parser.add_argument('--basis',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--start-row',type=int,default=0)
    parser.add_argument('--stop-row',type=int,default=0);parser.add_argument('--batch',type=int,default=2)
    parser.add_argument('--seconds',type=int,default=900);parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--points-override',type=int,default=0)
    run(**vars(parser.parse_args()))
