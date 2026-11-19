"""Materialize an isolated p=1000037 package from hash-pinned p=1000033 inputs.

The manifest contains science inputs and entry points, not output, binaries or old controls.
Usage: python prepare.py BASE EMPTY_DIRECTORY . The two trees must be separate.
"""
import argparse, hashlib, json, math, shutil
from pathlib import Path

P=1000037
HERE=Path(__file__).resolve().parent
BASE_MANIFEST='3cd48eea41430f6cd08449f651f729db8125417ce04b1d6b8dbd6f4b77867149'
FREE='a7421e462000b93c18376fbe5eda38c69e1cfa2b8642a7e5034e1824309b0e2d'
S='runtime/work/sp3-native.nosync/'
sha=lambda b: hashlib.sha256(b).hexdigest()

def sub(s,old,new,count=1):
    if s.count(old)!=count:raise ValueError((old,s.count(old),count))
    return s.replace(old,new)

def make(base,dest):
    assert __debug__ and all(P%i for i in range(2,math.isqrt(P)+1))
    assert base.is_dir() and base!=dest and base not in dest.parents and dest not in base.parents
    dest.mkdir(parents=True,exist_ok=True);assert not list(dest.iterdir()),'Requires empty destination'
    original=(base/'SOURCE_MANIFEST.json').read_bytes();assert sha(original)==BASE_MANIFEST
    original=json.loads(original); names=[n for n in original if n.split('/')[0] in ('runtime','first_basis','second_basis','quantum_source')];assert len(names)==447
    names+=['distributed_maps.py','solve_compact.py','build_native.py','compact_su3_maps_v1.py','validate_compact_su3_maps_v1.py','tasks.json']; digests={};edits={}
    for n in names:
        b=(base/n).read_bytes();assert sha(b)==original[n],n
        s=None
        if n in ('build_native.py',S+'native_su_projection_v1.py',S+'validate_su3_second_correction_map_v1.py',S+'validate_su3_second_quantum_map_v1.py'):
            s=sub(b.decode(),'(1000003,1000033)','(1000003,1000033,1000037)')
        elif n==S+'su_full_projection_native_v1.cpp':s=sub(b.decode(),'p!=1000003 && p!=1000033','p!=1000003 && p!=1000033 && p!=1000037')
        elif n in (S+'su3_correction_map_v2.py',S+'su3_second_correction_map_v1.py',S+'su3_coupled_solve_v1.py'):
            s=sub(b.decode(),'1000033','1000037')
        elif n==S+'validate_su3_second_coupled_solve_v1.py':
            s=sub(b.decode(),'checks=checks,equations=30,variables=16','prime=PRIME,checks=checks,equations=30,variables=16')
        elif n=='distributed_maps.py':
            s=sub(b.decode(),'prime=1000033','prime=1000037');s=sub(s,"assert json.loads((O/'preflight/result.json').read_text())['status']=='passed_independent_prime_runtime_maps_coupled_controls_and_Q1_json_sign'", "pre=json.loads((O/'preflight/result.json').read_text());assert pre['status']=='passed_p3' and pre['prime']==1000037 and pre['manifest_sha256']==hashlib.sha256((R/'SOURCE_MANIFEST.json').read_bytes()).hexdigest()")
            s=sub(s,'WALL=3300','WALL=1980')
        elif n=='solve_compact.py':
            s=sub(b.decode(),'kind,1000033','kind,1000037');s=sub(s,"for index in range(24):assert json.loads((O/('shard%d'%index)/'result.json').read_text())['status']=='complete_independent_prime_shard'","for index in range(24):\n    receipt=json.loads((O/('shard%d'%index)/'result.json').read_text());assert (receipt['prime'],receipt['shard'],receipt['status'])==(1000037,index,'complete_independent_prime_shard')")
            s+='''\nimport numpy as np\nrecord=json.loads((O/'second_solve_compact/result.json').read_text())\nif record['status']=='inconsistent_modular_second_correction_system':\n    raise RuntimeError('Inconsistent modular system: preserve stages and audit the obstruction')\nassert record['status']=='complete_modular_second_solution_with_all_first_parameters' and record['prime']==1000037\nfree=O/'second_solve_compact/free_variables.npy'; x=O/'second_solve_compact/particular.npy'\nassert hashlib.sha256(free.read_bytes()).hexdigest()==record['free_variables_sha256']==%r\nassert hashlib.sha256(x.read_bytes()).hexdigest()==record['particular_sha256']\nf=np.load(free,allow_pickle=False);v=np.load(x,allow_pickle=False)\nassert f.shape==(27492,) and v.shape==(77639,) and not np.any(v[f]) and np.all(v<1000037)\nassert record['equation_rank']==50147 and len(record['residuals'])==50 and all(r['nonzero']==0 for r in record['residuals'])\nassert [sum(r['coordinates'] for r in record['residuals'] if r['equation']==i) for i in (1,2)]==[34901,65406]\n'''%FREE
        if s is not None:b=s.encode();edits[n]=[original[n],sha(b)]
        dst=dest/n;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(b);digests[n]=sha(b)
    for n in ('preflight.py','postsolve.py'):
        b=(HERE/n).read_bytes();(dest/n).write_bytes(b);digests[n]=sha(b)
    tasks=json.loads((dest/'tasks.json').read_text());assert len(tasks)==24 and sum(map(len,tasks))==2197
    for kind,required in [('first_Q0',28311),('second_Q0',49328),('second_Q1',28311)]:
        intervals={}
        for shard in tasks:
            for item in shard:
                if item['kind']==kind:intervals.setdefault(item['index'],[]).append((item['start'],item['stop']))
        total=0
        for group in intervals.values():
            stop=0
            for a,b in sorted(group):assert a==stop<b;stop=b
            total+=stop
        assert total==required,(kind,total)
    for kind,count in [('first_rhs',13),('second_rhs',37)]:assert sorted(t['index'] for shard in tasks for t in shard if t['kind']==kind)==list(range(count))
    blob=(json.dumps(digests,sort_keys=True,indent=2)+'\n').encode();(dest/'SOURCE_MANIFEST.json').write_bytes(blob)
    return dict(prime=P,previous_manifest=BASE_MANIFEST,production_manifest=sha(blob),files=len(digests),prior_science_files_verified=447,modified_inputs=edits,shards=24,items=2197,ordered_free_coordinate_file_sha256=FREE)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('base',type=Path);p.add_argument('destination',type=Path);args=p.parse_args();print(json.dumps(make(args.base.resolve(),args.destination.resolve()),indent=2))
