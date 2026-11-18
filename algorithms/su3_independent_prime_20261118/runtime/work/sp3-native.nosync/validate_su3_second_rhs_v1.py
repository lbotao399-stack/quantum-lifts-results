"""All second RHS profiles: exact sign, chunking and independent matrices."""
from pathlib import Path
from fractions import Fraction
import argparse,json,tempfile,time
import su3_second_rhs_v1 as rhs
from sp_full_projected_basis import evaluate as dense
import numpy as np


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'))
    args=parser.parse_args();started=time.monotonic();expr=rhs.read_source(args.source)
    plan=json.loads((rhs.maps.HERE/'su3_second_correction_profiles_v1.json').read_text());by={}
    for word,c in expr.items():by.setdefault(rhs.profile(tuple(L for trace in word for L in trace)),{})[word]=c
    profiles={tuple((tuple(key),n) for key,n in row['profile']):row['index'] for row in plan['targets']}
    assert set(by)<=set(profiles)
    checks=[]
    with tempfile.TemporaryDirectory(dir=rhs.maps.HERE,prefix='second-rhs-check-') as name:
        root=Path(name)
        for row in plan['targets']:
            index=row['index'];part=by.get(tuple((tuple(key),n) for key,n in row['profile']),{})
            rhs.run(index,args.source,root/str(index),chunk=127,workers=2,points_override=2)
            frame,_=rhs.maps.target_frame(plan,index,points_override=2)
            value=dense([{w:-Fraction(c).numerator*pow(Fraction(c).denominator,-1,rhs.maps.PRIME)%rhs.maps.PRIME
                for w,c in part.items()}],frame)[0]
            assert np.array_equal(value,np.load(root/str(index)/'rhs.npy'))
            checks.append(dict(target=index,terms=len(part),coordinates=len(value),nonzero=int(np.count_nonzero(value))))
    result=dict(status='passed_all_second_rhs_profiles_negative_sign_chunks_and_dense_matrices',checks=checks,
        all_input_terms=sum(row['terms'] for row in checks),coordinates=sum(row['coordinates'] for row in checks),
        sources={name:rhs.maps.sha(rhs.maps.HERE/name) for name in
            ('su3_second_rhs_v1.py','validate_su3_second_rhs_v1.py','su3_second_correction_map_v1.py')},
        seconds=time.monotonic()-started)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key!='checks'}),flush=True)


if __name__=='__main__':main()
