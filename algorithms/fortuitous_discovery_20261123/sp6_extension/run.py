#!/usr/bin/env python3
"""Portable driver; install this directory inside the published base package."""
import argparse,importlib.util,json,os,subprocess,sys
from pathlib import Path

CODE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('part',choices=['build','small','solve','finalize','verify','reproduce','benchmark','next','recover-sp5','reuse'])
    ap.add_argument('--base',type=Path,default=CODE.parent)
    ap.add_argument('--workdir',type=Path,default=CODE/'work')
    ap.add_argument('--seconds',type=int,default=600)
    ap.add_argument('--max-sources',type=int,default=30000)
    args=ap.parse_args();base=args.base.resolve();work=args.workdir.resolve();work.mkdir(parents=True,exist_ok=True)
    assert (base/'replay.py').exists(),'Pass --base pointing to the complete published package.'
    spec=importlib.util.spec_from_file_location('base_package',base/'replay.py')
    baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)
    baseline.environment();baseline.hydrate(5)
    os.environ['FORTUITY_BASE_PACKAGE']=str(base)
    os.environ['FORTUITY_RUN_DIR']=str(work)
    env=dict(os.environ)
    env['PYTHONPATH']=str(CODE)+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    def execute(part):
        mapping={
          'build':('build_rank.py',['--rank','6','--seconds',str(args.seconds)]),
          'small':('small_complex.py',['--rank','6']),
          'solve':('support_solver.py',['--defect',str(work/'small_sp6/candidate_defect.pkl'),
                   '--output',str(work/'support_sp6'),'--max-sources',str(args.max_sources),
                   '--seconds',str(args.seconds)]),
          'finalize':('finalize_witness.py',['--rank','6']),
          'verify':('verify.py',[]),
          'benchmark':('benchmark_sp6.py',[]),
          'next':('next_rank_relations.py',[]),
          'recover-sp5':('recover_sp5.py',[]),
          'reuse':('discovery_scripts/spectator_reuse.py',[]),
        }
        script,extra=mapping[part]
        subprocess.run([sys.executable,str(CODE/script),*extra],env=env,check=True)
    if args.part=='reproduce':
        while not (work/'frames_sp6/construction.json').exists():execute('build')
        execute('small')
        while not (work/'support_sp6/result.json').exists():execute('solve')
        result=json.loads((work/'support_sp6/result.json').read_text())
        assert result['passed'],'A restricted inverse solve failed; inspect its recorded reason.'
        execute('finalize');execute('verify');execute('next')
    elif args.part=='verify':
        while not (work/'frames_sp6/construction.json').exists():execute('build')
        execute('verify')
    else:execute(args.part)

if __name__=='__main__':main()
