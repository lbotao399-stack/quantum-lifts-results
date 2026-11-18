"""Rebuild every saved -Q2 coordinate using the independent Python evaluator."""
import sp_environment
from pathlib import Path
from fractions import Fraction
from collections import Counter
import argparse,gzip,hashlib,json,time
import numpy as np
from fortuity.full.evalg import Frame

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def letters(raw):return [(x[0],x[1],tuple(x[2])) for x in raw]
def profile(word):return tuple(sorted(Counter((len(x[2]),x[0]+x[1]) for trace in word for x in trace).items()))

def run(evidence,output=None):
    started=time.monotonic();manifest=json.loads((evidence/'SHA256SUMS.json').read_text())
    for name,digest in manifest.items():
        part=Path(name);assert not part.is_absolute() and '..' not in part.parts
        assert sha(evidence/part)==digest,name
    plan=json.loads((evidence/'su3_second_correction_profiles_v1.json').read_text())
    source=evidence/'quantum_source/Q2.json.gz'
    declared=json.loads((source.parent/'result.json').read_text());assert sha(source)==declared['files'][source.name]
    with gzip.open(source,'rt') as stream:raw=json.load(stream)['terms']
    expression={tuple(tuple(letters(trace)) for trace in word):Fraction(c) for c,word in raw}
    assert len(expression)==len(raw)
    grouped={}
    for word,c in expression.items():grouped.setdefault(profile(word),{})[word]=c
    expected={tuple((tuple(key),n) for key,n in target['profile']) for target in plan['targets']}
    assert set(grouped)<=expected
    checks=[]
    for target in plan['targets']:
        index=target['index'];folder=evidence/'rhs'/('profile%02d'%index)
        declared=json.loads((folder/'result.json').read_text());input_record=json.loads((folder/'manifest.json').read_text())
        assert declared['status']=='complete_negative_Q2_evaluations' and declared['points_override']==0
        assert declared['target']==index and input_record['target']==target
        assert sha(folder/'rhs.npy')==declared['sha256']
        desc=input_record['frame'];prime=desc['prime'];assert prime==1000003 and desc['N']==3 and desc['preserved_target_basis_frame']
        frame=Frame(desc['N'],desc['points'],desc['k'],desc['degrees'],letters(desc['alphabet']),seed=desc['seed'],p=prime)
        if desc['additional_letters']:
            extra=Frame(desc['N'],desc['points'],desc['k'],desc['degrees'],letters(desc['additional_letters']),seed=desc['additional_seed'],p=prime)
            assert not set(frame.L).intersection(extra.L);frame.L.update(extra.L)
        assert frame.ncols==desc['columns']==declared['columns']
        part=grouped.get(tuple((tuple(key),n) for key,n in target['profile']),{})
        assert len(part)==input_record['input_terms']==declared['input_terms']
        words=sorted(part);value=np.zeros(frame.ncols,dtype=np.int64)
        for start in range(0,len(words),127):
            selected=words[start:start+127]
            coefficients=np.asarray([-part[w].numerator*pow(part[w].denominator,-1,prime)%prime for w in selected],dtype=np.int64)
            assert len(selected)*(prime-1)**2<2**63
            values=frame.eval_words(selected,dtype=np.int64)
            value=(value+coefficients@values)%prime
        saved=np.load(folder/'rhs.npy');assert np.array_equal(value,saved),index
        assert int(np.count_nonzero(value))==declared['nonzero_entries']
        checks.append(dict(target=index,terms=len(words),coordinates=frame.ncols,nonzero=int(np.count_nonzero(value))))
        print(json.dumps(checks[-1]),flush=True)
    assert sum(x['terms'] for x in checks)==len(expression)
    result=dict(status='passed_all_production_negative_Q2_coordinates_independent_Python_matrices',
        payload_hashes=len(manifest),profiles=len(checks),source_terms=len(expression),
        coordinates=sum(x['coordinates'] for x in checks),nonzero_coordinates=sum(x['nonzero'] for x in checks),
        checks=checks,seconds=time.monotonic()-started,source_json_sha256=sha(source),
        replay_driver_sha256=sha(Path(__file__)),
        scope='Complete replay of the supplied exact Q2 image in all production target frames. This does not solve the coupled equation or certify a quantum lift.')
    output=output or evidence/'replay.json';output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'}),flush=True);return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--evidence',type=Path,required=True);parser.add_argument('--output',type=Path)
    run(**vars(parser.parse_args()))
