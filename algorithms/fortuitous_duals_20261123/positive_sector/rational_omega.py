#!/usr/bin/env python3
"""Export the original rational omega and replay both integer closure frames."""
from probe_stable import HERE,m
from fractions import Fraction
from collections import Counter
import gzip,json,pickle,math,time,hashlib
from trace_descent import evaluate

def main():
    start=time.monotonic()
    path=HERE/'inputs/algorithms/dual_singlets_20261118/sp3_omega'
    load=lambda name:json.loads(gzip.decompress((path/name).read_bytes()))
    page=load('target.json.gz');record=load('cycle.json.gz');replay=load('cycle_replay.json.gz')
    _,centre,_,_=pickle.loads((HERE/'stable_rows.pkl').read_bytes())
    basis=page['centre_frames']['4']['minor_rows']
    coeffs=[(basis[i],Fraction(c)) for i,c in record['cycle_basis_coefficients']]
    denominator=math.lcm(*(c.denominator for _,c in coeffs))
    numerator=Counter()
    for i,c in coeffs:
        for w,n in centre[i].items():numerator[w]+=int(c*denominator)*n
    numerator={w:n for w,n in numerator.items() if n}
    variation=m.q(numerator)
    assert denominator==13824 and len(numerator)==4264 and len(variation)==4502
    checks=[]
    for label,frames in [('original',record['signed_frames']),('independent',replay['signed_frames'])]:
        frames={int(k):v for k,v in frames.items()}
        for frame in frames.values():
            frame['alphabet']=[(a,b,tuple(c)) for a,b,c in frame['alphabet']]
        values=evaluate([variation],frames,6,'Sp')[0]
        assert not any(values)
        checks.append({'label':label,'exact_integer_coordinates':len(values),'nonzero':0,'frame_dimensions':{str(k):v['dimension'] for k,v in frames.items()},'frame_minor_determinants_mod_prime':{str(k):v['determinant_mod_prime'] for k,v in frames.items()}})
        print(label,'exact zero',len(values),'seconds',time.monotonic()-start,flush=True)
    output={'weight':[3,3,5,0,0],'length':9,'coefficient_convention':'Original rational cycle.json.gz coefficients, no cochain rescaling','terms':[[str(Fraction(n,denominator)),m.encode(w)] for w,n in sorted(numerator.items())]}
    raw=(json.dumps(output,separators=(',',':'))+'\n').encode()
    (HERE/'omega_rational.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    result={'passed':True,'basis_coefficients':len(coeffs),'cochain_terms':len(numerator),'formal_derivative_terms':len(variation),'intermediate_common_denominator':denominator,'checks':checks,'omega_json_sha256':hashlib.sha256(raw).hexdigest(),'input_cycle_sha256':hashlib.sha256((path/'cycle.json.gz').read_bytes()).hexdigest(),'seconds':time.monotonic()-start,'scope':'Exact integer evaluation of the rational cycle variation on both complete archived frames; completeness uses their pinned rational character dimension and nonzero minor certificates.'}
    (HERE/'rational_omega_result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
