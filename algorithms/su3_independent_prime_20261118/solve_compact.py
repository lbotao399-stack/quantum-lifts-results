"""Pack completed maps and solve both equations with all first-order freedom."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

HERE=Path(__file__).resolve().parent
R=HERE;S=R/'runtime/work/sp3-native.nosync';O=R/'output'
assert os.environ.get('SLURM_JOB_ID');started=time.monotonic()
for name,digest in json.loads((HERE/'SOURCE_MANIFEST.json').read_text()).items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
for index in range(24):assert json.loads((O/('shard%d'%index)/'result.json').read_text())['status']=='complete_independent_prime_shard'
subprocess.run([sys.executable,str(HERE/'validate_compact_su3_maps_v1.py'),'--runtime',str(S),
    '--output',str(HERE/'controls.json')],check=True,timeout=120)
sys.path.insert(0,str(S))
import sp_environment
import compact_su3_maps_v1 as compact
import su3_coupled_solve_v1 as first
import su3_second_coupled_solve_v1 as second
from flint import ctx
ctx.threads=64
C=O/'compacted'
specs=[('first',R/'first_basis','su3_correction_profiles_v1.json','first_Q0'),
       ('second_Q0',R/'second_basis/shard0','su3_second_correction_profiles_v1.json','second_Q0'),
       ('second_Q1',R/'first_basis','su3_second_quantum_profiles_v1.json','second_Q1')]
records=[]
for name,basis,plan,kind in specs:
    result=compact.run(O/name,basis,C/name,json.loads((S/plan).read_text()),kind,1000033,workers=8)
    records.append({k:v for k,v in result.items() if k!='profiles_records'})
    print(json.dumps(records[-1]),flush=True)
for name,target in [('basis',R/'first_basis'),('rhs',O/'first/rhs')]:
    link=C/'first'/name
    if not link.exists():link.symlink_to(target,target_is_directory=True)
compact.save(C/'result.json',dict(status='complete_all_map_compaction',results=records,seconds=time.monotonic()-started))
first.run(C/'first',O/'first_solve_compact',seconds=max(1,int(6000-(time.monotonic()-started))))
record=json.loads((O/'first_solve_compact/result.json').read_text())
if record['status']!='complete_modular_evaluation_solution_with_all_free_parameters':
    print(json.dumps(record));sys.exit(2)
second.run(C/'first',O/'first_solve_compact',C/'second_Q0',C/'second_Q1',O/'second_rhs',
    [R/'second_basis'],O/'second_solve_compact',seconds=max(1,int(20000-(time.monotonic()-started))),threads=64)
