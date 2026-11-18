"""Check preserved frames and every structural block for selected profiles."""
from fractions import Fraction
from pathlib import Path
from math import ceil,comb
import argparse,json,random,time
import numpy as np
import su3_second_correction_map_v1 as maps
from sp_full_profiles import profile
from sp_full_projected_basis import evaluate as dense
from native_su_parallel_v1 import evaluate_profile
from su3_stream_rank_v3 import random_word


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args();started=time.monotonic()
    plan=json.loads((maps.HERE/'su3_second_correction_profiles_v1.json').read_text())
    mu=tuple(plan['mu']);length=plan['source_length'];all_contents=maps.contents(mu,length)
    indices=sorted({0,max(plan['sources'],key=lambda row:row['dimension'])['index'],
        max(plan['sources'],key=lambda row:max(key[1] for key,_ in row['profile']))['index'],
        max(plan['sources'],key=lambda row:row['fermions'])['index']})
    records=[]
    for index in indices:
        source=plan['sources'][index];selected=tuple((tuple(key),n) for key,n in source['profile'])
        choices=[content for content in all_contents if profile(content)==selected]
        rng=random.Random(2026111881+index)
        for attempt in range(2000):
            word=random_word(rng.choice(choices),rng,6)
            if word is None:continue
            initial=maps.algebra.parity_project({word:1},-1)
            if maps.project(initial,maps.sector(mu,length),1000003):break
        else:raise RuntimeError('No control word')
        rational=maps.algebra.projection.project(initial,maps.sector(mu,length))
        derivative=maps.algebra.differential(rational)
        assert not maps.algebra.differential(derivative)
        by_profile={}
        for w,c in derivative.items():
            key=profile(tuple(letter for trace in w for letter in trace));by_profile.setdefault(key,{})[w]=c
        allowed={tuple((tuple(key),n) for key,n in plan['targets'][j]['profile']) for j in source['targets']}
        assert set(by_profile)<=allowed
        for prime in (1000003,1000033):
            expression={w:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime for w,c in rational.items()}
            for j in source['targets']:
                frame,description=maps.target_frame(plan,j,prime,points_override=3)
                original=maps.Frame(3,3,description['k'],description['degrees'],description['alphabet'],seed=description['seed'],p=prime)
                assert all(np.array_equal(frame.L[L][1],original.L[L][1]) for L in original.L)
                target_profile=tuple((tuple(key),n) for key,n in plan['targets'][j]['profile'])
                explicit={w:Fraction(c).numerator*pow(Fraction(c).denominator,-1,prime)%prime
                    for w,c in by_profile.get(target_profile,{}).items()}
                expected=dense([explicit],frame)
                for workers in (1,2):
                    actual=evaluate_profile([expression],frame,source['profile'],plan['targets'][j]['profile'],workers)
                    assert np.array_equal(actual,expected),(index,j,prime,workers)
                records.append(dict(source=index,target=j,prime=prime,coordinates=expected.size,
                    nonzero=int(np.count_nonzero(expected)),original_target_matrices_preserved=True,
                    additional_letters=len(description['additional_letters'])))
    assert any(row['nonzero'] for row in records)
    result=dict(status='passed_second_Q0_profile_blocks_rational_dense_and_frame_preservation',
        checks=records,source_profiles=indices,seconds=time.monotonic()-started,
        sources={name:maps.sha(maps.HERE/name) for name in
            ('su3_second_correction_map_v1.py','validate_su3_second_correction_map_v1.py','native_su_parallel_v1.py')})
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],checks=len(records),coordinates=sum(row['coordinates'] for row in records),seconds=result['seconds'])),flush=True)


if __name__=='__main__':main()
