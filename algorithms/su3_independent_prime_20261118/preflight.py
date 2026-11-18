from pathlib import Path
import gzip,hashlib,json,pickle,subprocess,sys,tempfile,time

R=Path(__file__).resolve().parent;S=R/'runtime/work/sp3-native.nosync'
started=time.monotonic()
for name,digest in json.loads((R/'SOURCE_MANIFEST.json').read_text()).items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
out=R/'output/preflight';out.mkdir(parents=True,exist_ok=True)
commands=[['build_native.py'],[str(S/'validate_su3_second_correction_map_v1.py')],
    [str(S/'validate_su3_second_quantum_map_v1.py'),'--basis',str(R/'first_basis'),'--output',str(out/'Q1_controls.json')],
    [str(S/'validate_su3_second_coupled_solve_v1.py')]]
for index,args in enumerate(commands):
    with (out/('control%d.log'%index)).open('w') as log:
        subprocess.run([sys.executable,*args],cwd=R,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
sys.path.insert(0,str(S))
import su3_correction_map_v2 as first
import su3_correction_rhs_v1 as rhs
import su3_second_correction_map_v1 as second
import su3_second_quantum_map_v1 as quantum
import su3_coupled_solve_v1 as solve
import numpy as np
from fractions import Fraction
from sp_full_profiles import profile
from sp_full_projected_basis import evaluate as dense
assert first.PRIME==rhs.PRIME==second.PRIME==quantum.PRIME==solve.PRIME==1000033
plan=json.loads((S/'su3_correction_profiles_v1.json').read_text())
with gzip.open(R/'quantum_source/Q1.json.gz','rt') as stream:raw=json.load(stream)['terms'][:2]
expression={tuple(tuple((letter[0],letter[1],tuple(letter[2])) for letter in trace) for trace in word):Fraction(c) for c,word in raw}
key=profile(tuple(letter for trace in next(iter(expression)) for letter in trace))
target=next(t for t in plan['targets'] if tuple((tuple(k),n) for k,n in t['profile'])==key)
with tempfile.TemporaryDirectory(dir=out) as temp:
    folder=Path(temp);source=folder/'Q1.json.gz'
    with gzip.open(source,'wt') as stream:json.dump(dict(terms=raw),stream)
    (folder/'result.json').write_text(json.dumps(dict(files={source.name:first.sha(source)})))
    rhs.run(target['index'],source,folder/'evaluated',workers=2)
    frame,description=first.target_frame(target)
    expected=dense([{w:-c.numerator*pow(c.denominator,-1,first.PRIME)%first.PRIME for w,c in expression.items() if profile(tuple(letter for trace in w for letter in trace))==key}],frame)[0]
    assert np.array_equal(expected,np.load(folder/'evaluated/rhs.npy'))
result=dict(status='passed_independent_prime_runtime_maps_coupled_controls_and_Q1_json_sign',prime=1000033,
    first_rhs_dense_coordinates=len(expected),seconds=time.monotonic()-started,
    manifest_sha256=hashlib.sha256((R/'SOURCE_MANIFEST.json').read_bytes()).hexdigest())
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
