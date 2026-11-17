"""All first-profile Q1 supports, rational reduction and independent matrices."""
from fractions import Fraction
from pathlib import Path
import argparse,json,pickle,tempfile,time
import numpy as np
import su3_second_quantum_map_v1 as driver
from sp_full_profiles import profile
from sp_full_projected_basis import evaluate as dense
from native_su_parallel_v1 import evaluate


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--basis',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args();started=time.monotonic();maps=driver.maps
    plan=json.loads((driver.HERE/'su3_second_quantum_profiles_v1.json').read_text())
    second=json.loads((driver.HERE/'su3_second_correction_profiles_v1.json').read_text())
    records=[];reference_controls=[]
    for source in plan['sources']:
        index=source['index'];folder=args.basis/('profile%02d'%index)
        declared=json.loads((folder/'result.json').read_text())
        assert maps.sha(folder/'independent_seeds.pkl')==declared['seeds_sha256']
        seeds=pickle.loads((folder/'independent_seeds.pkl').read_bytes());position=len(seeds)//2;word=seeds[position]
        initial=maps.algebra.parity_project({word:1},-1)
        rational=maps.algebra.projection.project(initial,maps.sector(tuple(plan['source_mu']),9))
        exact=driver.quantum.act(rational,1,3)
        if index in (0,6,16):
            part=dict(sorted(rational.items())[:3])
            assert driver.quantum.act(part,1,3)==driver.quantum.reference(part,1,3)
            reference_controls.append(dict(profile=index,input_words=len(part)))
        allowed={tuple((tuple(key),n) for key,n in plan['targets'][j]['profile']):j for j in source['targets']}
        for prime in (1000003,1000033):
            _,actual=driver.quantum_image(word,tuple(plan['source_mu']),prime)
            expected={w:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime for w,c in exact.items()}
            expected={w:c for w,c in expected.items() if c}
            assert actual==expected,(index,prime)
            parts={j:{} for j in source['targets']}
            for w,c in expected.items():
                key=profile(tuple(letter for trace in w for letter in trace))
                assert key in allowed,(index,key)
                parts[allowed[key]][w]=c
            for j in source['targets']:
                frame,_=maps.target_frame(second,j,prime,points_override=2)
                reference=dense([parts[j]],frame)
                assert np.array_equal(reference,evaluate([parts[j]],frame,2)),(index,j,prime)
                records.append(dict(source=index,source_row=position,target=j,prime=prime,
                    quantum_terms=len(parts[j]),coordinates=reference.size,nonzero=int(np.count_nonzero(reference))))
        print(json.dumps(dict(stage='profile_checked',source=index,rational_terms=len(rational),quantum_terms=len(exact))),flush=True)
    # Exercise nonzero global row offsets and the on-disk production assembly.
    with tempfile.TemporaryDirectory(prefix='second-Q1-control-',dir=driver.HERE) as name:
        output=Path(name);result=driver.run(0,args.basis/'profile00',output,start_row=1,stop_row=3,batch=1,workers=2,points_override=2)
        seeds=pickle.loads((args.basis/'profile00/independent_seeds.pkl').read_bytes())
        for row in result['records']:
            frame,_=maps.target_frame(second,row['target'],points_override=2)
            selected=tuple((tuple(key),n) for key,n in plan['targets'][row['target']]['profile'])
            expression=[]
            for word in seeds[row['start']:row['stop']]:
                _,image=driver.quantum_image(word,tuple(plan['source_mu']))
                expression.append({w:c for w,c in image.items() if profile(tuple(letter for trace in w for letter in trace))==selected})
            expected=dense(expression,frame)
            path=output/('t%02d_%05d_%05d.npy'%(row['target'],row['start'],row['stop']))
            assert np.array_equal(np.load(path),expected)
    result=dict(status='passed_all_Q1_profile_supports_rational_reductions_dense_matrices_and_global_rows',
        checks=records,reference_controls=reference_controls,structural_edges=plan['structural_edges'],
        coordinates=sum(row['coordinates'] for row in records),source_profiles=len(plan['sources']),
        sources={name:maps.sha(driver.HERE/name) for name in ('su3_second_quantum_profiles_v1.py',
            'su3_second_quantum_profiles_v1.json','su3_second_quantum_map_v1.py','validate_su3_second_quantum_map_v1.py')},
        seconds=time.monotonic()-started)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ('checks','sources')}),flush=True)


if __name__=='__main__':main()
