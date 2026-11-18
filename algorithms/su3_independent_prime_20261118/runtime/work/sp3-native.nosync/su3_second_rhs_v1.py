"""Pinned rational -Q2 x in the preserved nonabelian second-target frames."""
from fractions import Fraction
from pathlib import Path
import argparse,gzip,json,pickle,time
import numpy as np
import su3_second_correction_map_v1 as maps
from sp_full_profiles import profile
from native_su_parallel_v1 import evaluate


def read_source(source):
    declared=json.loads((source.parent/'result.json').read_text())
    assert source.name in ('Q2.pkl','Q2.json.gz')
    assert maps.sha(source)==declared['files'][source.name]
    if source.name=='Q2.pkl':return pickle.loads(source.read_bytes())['expression']
    with gzip.open(source,'rt') as stream:terms=json.load(stream)['terms']
    return {tuple(tuple((letter[0],letter[1],tuple(letter[2])) for letter in trace) for trace in word):Fraction(c)
            for c,word in terms}


def run(index,source,output,chunk=512,workers=8,points_override=0):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=True)
    plan=json.loads((maps.HERE/'su3_second_correction_profiles_v1.json').read_text())
    target=plan['targets'][index];selected=tuple((tuple(key),n) for key,n in target['profile'])
    expression=read_source(source)
    terms=[(word,-Fraction(c).numerator*pow(Fraction(c).denominator,-1,maps.PRIME)%maps.PRIME)
        for word,c in sorted(expression.items()) if profile(tuple(letter for trace in word for letter in trace))==selected]
    frame,description=maps.target_frame(plan,index,points_override=points_override)
    manifest=dict(target=target,frame=description,source_sha256=maps.sha(source),input_terms=len(terms),
        all_input_terms=len(expression),chunk=chunk,points_override=points_override,
        driver_sha256=maps.sha(__file__),frame_driver_sha256=maps.sha(maps.__file__))
    manifest=json.loads(json.dumps(manifest));marker=output/'manifest.json'
    if marker.exists():assert json.loads(marker.read_text())==manifest
    else:maps.save(marker,manifest)
    total=np.zeros(frame.ncols,dtype=np.int64)
    for start in range(0,len(terms),chunk):
        path=output/('chunk%06d.npy'%start);receipt=path.with_suffix('.json')
        if receipt.exists():
            record=json.loads(receipt.read_text());assert maps.sha(path)==record['sha256'];value=np.load(path).astype(np.int64)
        else:
            before=time.monotonic();value=evaluate([dict(terms[start:start+chunk])],frame,workers)[0]
            temp=path.with_suffix('.tmp.npy');np.save(temp,value.astype(np.uint32));temp.replace(path)
            maps.save(receipt,dict(start=start,stop=min(start+chunk,len(terms)),seconds=time.monotonic()-before,sha256=maps.sha(path)))
        total=(total+value)%maps.PRIME
        maps.save(output/'progress.json',dict(target=index,completed_terms=min(start+chunk,len(terms)),required_terms=len(terms),seconds=time.monotonic()-started))
    path=output/'rhs.npy';np.save(path,total.astype(np.uint32))
    result=dict(status='complete_negative_Q2_evaluations',target=index,input_terms=len(terms),columns=frame.ncols,
        nonzero_entries=int(np.count_nonzero(total)),sha256=maps.sha(path),source_sha256=maps.sha(source),
        seconds=time.monotonic()-started,points_override=points_override,
        scope='The exact -Q2 source evaluated in the second-target frames. A second correction and the full Q1 homogeneous quotient remain unsolved.')
    maps.save(output/'result.json',result);print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--chunk',type=int,default=512);parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--points-override',type=int,default=0)
    run(**vars(parser.parse_args()))
