"""Twenty-four disjoint shards, each using four processes and eight native threads."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import fcntl,hashlib,json,os,subprocess,sys,time

SCHEDULE=Path(__file__).resolve().parent;R=SCHEDULE;S=R/'runtime/work/sp3-native.nosync';O=R/'output'
INDEX=int(os.environ['SLURM_ARRAY_TASK_ID']);assert 0<=INDEX<24
assert int(os.environ['SLURM_CPUS_PER_TASK'])==32
START=time.monotonic();WALL=3300
def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix('.tmp');temp.write_text(json.dumps(data,indent=2)+'\n');temp.replace(path)
for name,digest in json.loads((R/'SOURCE_MANIFEST.json').read_text()).items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==digest,name
assert json.loads((O/'preflight/result.json').read_text())['status']=='passed_independent_prime_runtime_maps_coupled_controls_and_Q1_json_sign'
for name,digest in json.loads((SCHEDULE/'SOURCE_MANIFEST.json').read_text()).items():
    assert hashlib.sha256((SCHEDULE/name).read_bytes()).hexdigest()==digest,name
tasks=json.loads((SCHEDULE/'tasks.json').read_text())[INDEX]
label='job'+os.environ['SLURM_JOB_ID']
def task(item):
    kind,index=item['kind'],item['index']
    if time.monotonic()-START>WALL-90:return dict(item,status='not_started_at_global_guard')
    limit=min(1800,max(1,int(WALL-45-(time.monotonic()-START))))
    if kind.endswith('rhs'):
        base=O/('first' if kind=='first_rhs' else 'second_rhs')
        output=base/'rhs'/('profile%02d'%index)
        name='su3_correction_rhs_v1.py' if kind=='first_rhs' else 'su3_second_rhs_v1.py'
        source=R/'quantum_source'/('Q1.json.gz' if kind=='first_rhs' else 'Q2.json.gz')
        args=['--index',str(index),'--source',str(source),'--output',str(output),'--workers','8']
    else:
        output=O/('first' if kind=='first_Q0' else kind)/'maps'/('profile%02d'%index)/('rows%05d_%05d'%(item['start'],item['stop']))
        if kind=='first_Q0':name='su3_correction_map_v2.py';basis=R/'first_basis'
        elif kind=='second_Q0':name='su3_second_correction_map_v1.py';basis=R/'second_basis/shard0'/('sources%02d'%index)
        else:name='su3_second_quantum_map_v1.py';basis=R/'first_basis'/('profile%02d'%index)
        args=['--index',str(index),'--basis',str(basis),'--output',str(output),'--start-row',str(item['start']),
            '--stop-row',str(item['stop']),'--batch','2' if kind=='second_Q1' else '8','--workers','8','--seconds',str(limit)]
    output.mkdir(parents=True,exist_ok=True)
    with (output/'execution.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        marker=output/'result.json';previous=json.loads(marker.read_text()) if marker.exists() else None
        if previous and previous['status'].startswith('complete_'):
            return dict(item,status='complete_previously',seconds=0)
        before=time.monotonic();peak=0;reason=None
        with (output/(label+'.log')).open('w') as log:
            proc=subprocess.Popen([sys.executable,'-u',str(S/name),*args],cwd=S,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
            while proc.poll() is None:
                try:
                    lines=(Path('/proc')/str(proc.pid)/'status').read_text().splitlines()
                    peak=max(peak,next((int(x.split()[1])*1024 for x in lines if x.startswith('VmRSS:')),0))
                except FileNotFoundError:pass
                if peak>12*1024**3:reason='memory_guard'
                if time.monotonic()-before>limit+30:reason='time_guard'
                if reason:
                    proc.terminate()
                    try:proc.wait(timeout=5)
                    except subprocess.TimeoutExpired:proc.kill();proc.wait()
                    break
                time.sleep(.5)
        result=json.loads(marker.read_text()) if marker.exists() else None
        record=dict(item,status=result['status'] if result else 'missing_result',exit_code=proc.returncode,
            stop_reason=reason,peak_RSS=peak,seconds=time.monotonic()-before)
        save(output/(label+'.receipt.json'),record);return record
finished=[];folder=O/('shard%d'%INDEX)
with ThreadPoolExecutor(max_workers=4) as pool:
    for future in as_completed([pool.submit(task,item) for item in tasks]):
        finished.append(future.result());save(folder/'progress.json',dict(finished=finished,required=len(tasks)))
complete=all(x['status'].startswith('complete_') and not x.get('exit_code') and not x.get('stop_reason') for x in finished)
record=dict(status='complete_independent_prime_shard' if complete else 'partial_independent_prime_shard',
    prime=1000033,job=os.environ['SLURM_JOB_ID'],node=os.environ.get('SLURMD_NODENAME'),shard=INDEX,
    cpus=32,processes=4,native_threads_per_process=8,seconds=time.monotonic()-START,finished=finished,
    d1_increment=0,new_generator_increment=0,d2_increment=0,d3_increment=0,np_increment=0)
save(folder/'result.json',record);print(json.dumps({k:v for k,v in record.items() if k!='finished'}),flush=True)
if not complete:sys.exit(2)
