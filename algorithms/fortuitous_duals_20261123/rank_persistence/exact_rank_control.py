"""Exact rational evaluations of the archived omega's stable derivative."""
from pathlib import Path
from fractions import Fraction
import json,time,hashlib,gzip
from rank_runtime import HERE,POSITIVE,m,load_omega
from terminal_interval_ranks import sample

def main():
    started=time.monotonic()
    omega=load_omega()
    delta=m.q(omega)
    assert len(delta)==4502
    words=sorted(delta)
    alphabet=sorted({l for w in words for tr in w for l in tr})
    records=[]
    for seed in (159771,159772,159773,159774):
        for size in (2,4,6,8):
            values=sample(words,5,size,1,seed,'Sp',alphabet)
            pairing=[sum(delta[w]*int(values[i,j]) for i,w in enumerate(words))
                     for j in range(values.shape[1])]
            row={'defining_dimension':size,'points':1,'seed':seed,
                 'exact_derivative_coordinates':[str(x) for x in pairing]}
            if size<=6:assert not any(pairing),row
            records.append(row);print(row,flush=True)
        if any(pairing):break
    assert any(pairing),'no nonzero larger-rank control in the tested samples'
    result={'passed':True,'source_commit':'9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9',
            'weight':[3,3,5,0,0],'cycle_trace_terms':len(omega),
            'stable_derivative_terms':len(delta),'alphabet':alphabet,'records':records,
            'integer_frame':'Pinned terminal_interval_ranks.sample; each matrix has integer entries and row L1 norm at most2; all final sums use Fraction.',
            'scope':'These are exact sample values and a nonzero chosen-lift witness at Sp4. Completeness of rank disappearance is supplied separately by the full1198-minor certificate, not by these samples.',
            'omega_json_sha256':hashlib.sha256(gzip.decompress((POSITIVE/'omega_rational.json.gz').read_bytes())).hexdigest(),
            'omega_gzip_sha256':hashlib.sha256((POSITIVE/'omega_rational.json.gz').read_bytes()).hexdigest(),
            'seconds':time.monotonic()-started}
    (HERE/'exact_rank_control.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
