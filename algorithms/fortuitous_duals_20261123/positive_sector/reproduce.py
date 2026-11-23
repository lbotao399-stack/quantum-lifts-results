#!/usr/bin/env python3
"""Run the complete Sp(3) finite-rank identity and fortuity calculation."""
from pathlib import Path
import argparse,subprocess,sys,os,json,time,hashlib

HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    subprocess.run([sys.executable,str(HERE/'fetch_inputs.py'),'--repo',str(args.repo.resolve())],check=True)
    manifest=json.loads((HERE/'inputs/INPUT_MANIFEST.json').read_text())
    for item in manifest:
        assert hashlib.sha256((HERE/'inputs'/item['path']).read_bytes()).hexdigest()==item['sha256']
    if args.prepare_only:
        print(json.dumps({'prepared':True,'verified_input_files':len(manifest)}));return
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1')
    env['N4_DUAL_SOURCE']=str(HERE/'inputs/algorithms/dual_singlets_20261118')
    env['N4_TRACE_RUNTIME']=str(HERE/'inputs/algorithms/su3_independent_prime_20261118/runtime/work/sp3-native.nosync')
    cache=HERE/'stable_rows.pkl'
    if cache.exists():cache.unlink()
    records=[]
    for name in ['probe_stable.py','build_native_frames.py','relation_image.py','rational_omega.py','certify_rational_relations.py','exact_relation_image.py']:
        start=time.monotonic()
        subprocess.run([sys.executable,str(HERE/name)],cwd=HERE,env=env,check=True)
        records.append({'script':name,'passed':True,'seconds':time.monotonic()-start})
    (HERE/'reproduction_run.json').write_text(json.dumps({'passed':True,'verified_inputs':len(manifest),'checks':records},indent=2)+'\n')

if __name__=='__main__':main()
