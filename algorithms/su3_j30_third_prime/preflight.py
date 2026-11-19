"""Build and exercise the production maps, RHS signs, packing and affine solver before admission."""
from pathlib import Path
from fractions import Fraction
import gzip,hashlib,json,os,pickle,subprocess,sys,tempfile,time
R=Path(__file__).resolve().parent; S=R/'runtime/work/sp3-native.nosync'; out=R/'output/preflight';out.mkdir(parents=True,exist_ok=True)
started=time.monotonic();P=1000037; sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest(); load=lambda p:json.loads(Path(p).read_bytes())
assert __debug__
for name,digest in load(R/'SOURCE_MANIFEST.json').items(): assert sha(R/name)==digest,name
commands=[('build',[R/'build_native.py']),('Q0',[S/'validate_su3_second_correction_map_v1.py','--output',out/'Q0.json']),('Q1',[S/'validate_su3_second_quantum_map_v1.py','--basis',R/'first_basis','--output',out/'Q1.json']),('coupled',[S/'validate_su3_second_coupled_solve_v1.py']),('pack',[R/'validate_compact_su3_maps_v1.py','--runtime',S,'--output',out/'pack.json']),('negative_Q2',[S/'validate_su3_second_rhs_v1.py','--source',R/'quantum_source/Q2.json.gz','--output',out/'negative_Q2.json'])]
for label,args in commands:
    with (out/(label+'.log')).open('w') as log:subprocess.run([sys.executable,*map(str,args)],cwd=R,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    print(label,'passed',flush=True)
sys.path.insert(0,str(S));import su3_correction_map_v2 as first,su3_correction_rhs_v1 as rhs,su3_second_correction_map_v1 as second,su3_second_quantum_map_v1 as quantum,su3_coupled_solve_v1 as solve
from sp_full_profiles import profile
from sp_full_projected_basis import evaluate as dense
import numpy as np
assert first.PRIME==rhs.PRIME==second.PRIME==quantum.PRIME==solve.PRIME==P
plan=load(S/'su3_correction_profiles_v1.json')
with gzip.open(R/'quantum_source/Q1.json.gz','rt') as f:raw=json.load(f)['terms'][:2]
expr={tuple(tuple((l[0],l[1],tuple(l[2])) for l in trace) for trace in word):Fraction(c) for c,word in raw}
key=profile(tuple(l for t in next(iter(expr)) for l in t)); target=next(t for t in plan['targets'] if tuple((tuple(k),n) for k,n in t['profile'])==key)
with tempfile.TemporaryDirectory(dir=out) as name:
    path=Path(name);file=path/'Q1.json.gz'
    with gzip.open(file,'wt') as f:json.dump(dict(terms=raw),f)
    (path/'result.json').write_text(json.dumps(dict(files={file.name:sha(file)})))
    rhs.run(target['index'],file,path/'rhs',workers=2)
    frame,_=first.target_frame(target);expected=dense([{w:-c.numerator*pow(c.denominator,-1,P)%P for w,c in expr.items() if profile(tuple(l for t in w for l in t))==key}],frame)[0]
    assert np.count_nonzero(expected)>0 and np.array_equal(expected,np.load(path/'rhs/rhs.npy'))
q0=[r for r in load(out/'Q0.json')['checks'] if r['prime']==P];q1=[r for r in load(out/'Q1.json')['checks'] if r['prime']==P];r2=load(out/'negative_Q2.json');pack=load(out/'pack.json');couple=load(S/'validate_su3_second_coupled_solve_v1.json');assert len({r['source'] for r in q1})==17 and max(r['nonzero'] for r in q0)>0 and max(r['nonzero'] for r in q1)>0 and len(r2['checks'])==37 and r2['all_input_terms']==23055 and sum(r['nonzero'] for r in r2['checks'])>0
assert pack['prime']==P and pack['matrix_coordinates_compared']==126 and couple['prime']==P
result=dict(status='passed_p3',prime=P,manifest_sha256=sha(R/'SOURCE_MANIFEST.json'),first_rhs_nonzero_coordinates=int(np.count_nonzero(expected)),first_rhs_coordinates=len(expected),second_rhs_exact_terms=r2['all_input_terms'],second_rhs_targets=len(r2['checks']),second_rhs_coordinates=r2['coordinates'],Q0_edges=len(q0),Q0_coordinates=sum(r['coordinates'] for r in q0),Q1_edges=len(q1),Q1_profiles=17,Q1_coordinates=sum(r['coordinates'] for r in q1),packed_coordinates=126,coupled_control=couple,seconds=time.monotonic()-started)
part=out/'result.tmp';part.write_text(json.dumps(result,indent=2)+'\n');part.replace(out/'result.json');print(json.dumps(result),flush=True)
