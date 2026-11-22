"""Rebuild the mixed SO(13) third pairing in a fresh temporary directory."""
import argparse
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def replay(repo, workers, lower):
    root = Path(__file__).resolve().parent
    repo = repo.resolve()
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    if 'PYTHONPATH' in env:
        env['PYTHONPATH'] = os.pathsep.join(str(Path(p).resolve()) for p in env['PYTHONPATH'].split(os.pathsep))
    with tempfile.TemporaryDirectory() as scratch:
        run_dir = Path(scratch)
        def run(*args):
            subprocess.run(list(map(str,args)),cwd=run_dir,env=env,check=True)
        def script(name,*args):
            run(sys.executable,root/name,'--repo',repo,*args)
        def same(name):
            assert (root/name).read_bytes() == (run_dir/name).read_bytes(), name
        for name in ('bmn_source_pair','bmn_pullback','bmn_forward_log'):
            run('c++','-std=c++17','-O3',root/(name+'.cpp'),'-o',run_dir/name)
        dual = root/'mixed13-dual.json'
        rows = root/'mixed13-rows.json'
        if lower:
            script('mixed_second.py','--colours',13,'--power',6,'--scalars',2,0,0,'--workers',workers,
                   '--out','mixed13-dual.json')
            same('mixed13-dual.json')
            script('mixed_adjoint_check.py','--dual',dual,'--workers',workers,'--out','mixed13-adjoints.json')
            same('mixed13-adjoints.json')
            script('mixed_relations.py','--dual',dual,'--workers',workers,'--out','mixed13-rows.json')
            same('mixed13-rows.json')
        seed = 'bmn13-third-seed.dat'
        script('bmn_third_dual.py','--dual',dual,'--rows',rows,'--out',seed)
        same(seed)
        tables = []
        parents = []
        for index,prime in enumerate((1000003,1000033),1):
            table = f'dual-{prime}.dat'
            parent = f'parents-{prime}.dat'
            run(run_dir/'bmn_pullback',prime,seed,table)
            run(run_dir/'bmn_source_pair',prime,table,root/'bmn13-parameters.dat',2,parent,'--all')
            tables.append(table)
            parents.append(parent)
        script('bmn_exact_dual.py','--seed',seed,'--modular',*tables,'--out','bmn13-exact-dual.dat')
        same('bmn13-exact-dual.dat')
        run(run_dir/'bmn_forward_log',parents[0],root/'bmn13-parameters.dat','active.dat')
        script('bmn_rational_pair.py','--parameters',root/'bmn13-parameters.dat','--active','active.dat',
               '--dual','bmn13-exact-dual.dat','--parents',*parents,'--workers',workers,
               '--out','bmn13-rational-pair.json')
        same('bmn13-rational-pair.json')
        script('bmn_substitution_check.py','--case',rows,'--out','bmn13-substitution-check.json')
        same('bmn13-substitution-check.json')
        print('verified exact pairing and substitution identities',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--rebuild-adjoints', action='store_true')
    a=p.parse_args()
    assert __debug__ and a.workers>=1
    replay(a.repo,a.workers,a.rebuild_adjoints)
