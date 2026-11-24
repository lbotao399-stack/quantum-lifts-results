#!/usr/bin/env python3
"""Execute the complete physical block certificate, with explicit cache scope."""
from pathlib import Path
import argparse,datetime,hashlib,importlib.metadata,json,os,subprocess,sys,time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PREFIX='algorithms/su3_second_correction_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor/fortuity/'

def save(name,value):
    temporary=HERE/(name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(HERE/name)

def source_manifest():
    base=json.loads((HERE/'declaration.json').read_text())['base_commit']
    sources=[]
    for name in ('__init__.py','full/__init__.py','full/letters.py','full/words.py','full/wordenum.py'):
        relative=PREFIX+name;path=ROOT/relative
        expected=subprocess.check_output(['git','show',f'{base}:{relative}'],cwd=ROOT)
        assert path.read_bytes()==expected,relative
        sources.append({'path':relative,'sha256':hashlib.sha256(expected).hexdigest()})
    save('input_manifest.json',{'public_repository':'https://github.com/lbotao399-stack/quantum-lifts-results',
          'base_commit':base,'all_imported_repository_sources':sources,
          'primary_sources':[{'url':'https://arxiv.org/pdf/1305.6314v2','equations':'2.7--2.11'},
                             {'url':'https://arxiv.org/pdf/2306.04673v2','equations':'A.1--A.5,D.4,D.9,3.13,E.21--E.30',
                              'normalization_audit':'The f primary factor in E.12/E.13 is not used; direct radial curvature contraction gives the value consistent with 3.13/E.30.'}],
          'python':sys.version,'packages':{name:importlib.metadata.version(name) for name in ('python-flint','numpy','sympy')}})

def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--full',action='store_true');mode.add_argument('--verify',action='store_true')
    args=parser.parse_args();source_manifest();started=time.monotonic()
    scope='full regeneration of all native frames and all Wick matrices' if args.full else 'exact certificate replay using stored Wick matrices'
    # A killed replay must never leave an earlier success receipt masquerading
    # as the status of the current files.
    save('execution_receipt.json',{'passed':False,'status':'running',
          'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'mode':scope,'steps':[]})
    if args.full:
        for pattern in ('frames_Y*.json.gz','full_gram_Y*.json.gz','finite_metric_Y*.json.gz',
                        'outgoing_gram_Y*.json.gz','harmonic_Y*.json.gz','physical_representatives.json.gz','physical_metric.json.gz'):
            for path in HERE.glob(pattern):path.unlink()
    environment=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    commands=[['build_frames.py']]+[['compute_metric.py','--length',str(n)] for n in (5,6,7,8)]
    commands +=[[name] for name in ('character_check.py','check_fortuity.py','independent_checks.py','export_block.py')]
    records=[]
    with (HERE/'execution.log').open('w') as output:
        for command in commands:
            step=time.monotonic();print('Running',*command,flush=True)
            result=subprocess.run([sys.executable,str(HERE/command[0]),*command[1:]],cwd=ROOT,env=environment,
                                  stdout=output,stderr=subprocess.STDOUT)
            output.flush()
            records.append({'command':command,'returncode':result.returncode,'seconds':time.monotonic()-step})
            save('execution_receipt.json',{'passed':False,'status':'running' if not result.returncode else 'failed',
                  'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'mode':scope,'steps':records,'seconds':time.monotonic()-started})
            if result.returncode:raise RuntimeError(command)
    generated=[]
    for path in sorted(HERE.iterdir()):
        if path.suffix in ('.json','.gz') and path.name not in ('execution_receipt.json','MANIFEST.json'):
            generated.append({'path':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    save('execution_receipt.json',{'passed':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'mode':scope,'status':'complete',
          'steps':records,'generated_file_hashes':generated,'seconds':time.monotonic()-started})
    print('Physical block certificate passed;',len(records),'steps;',time.monotonic()-started,'seconds.',flush=True)

if __name__=='__main__':main()
