"""Reconstruct three exact rank-lower witness matrices with the pinned engine."""
from pathlib import Path
from hashlib import sha256
import argparse, json, os, subprocess, tempfile

p=argparse.ArgumentParser()
p.add_argument('native_dir',type=Path)
a=p.parse_args()
root=a.native_dir.resolve()
f=root/'sha256.json'
assert sha256(f.read_bytes()).hexdigest()=='e9b4b27cc87fc4970ea43b65db417227695c90c9569a8d065f7c660bfec8f525'
for name,h in json.loads(f.read_bytes()).items():
    assert sha256((root/name).read_bytes()).hexdigest()==h,(name,h)
slots=[
    ((2,1,0,1,1,4,4),(12,4,4,4)),
    ((2,2,1,1,1,5,5),(101,46,49,6)),
    ((1,1,0,2,2,4,4),(28,8,15,5)),
]
env=dict(os.environ)
env.pop('A_NORHO',None)
env.update(PAR_THREADS='2',A_CHECK='1',A_VERB='0',A_SCALAR='1')
with tempfile.TemporaryDirectory(prefix='r007-small-rank-') as t:
    exe=Path(t)/'engine'
    subprocess.run(['c++','-O3','-std=c++17','-pthread',str(root/'engine.cpp'),'-o',str(exe)],check=True)
    for prime,seed in ((1000003,38221),(2147483647,38221),(1000003,38117)):
        for spec,expected in slots:
            r=json.loads(subprocess.check_output([str(exe),'chain','SU3',*map(str,(*spec,seed,prime))],env=env))
            assert r['sector']==list(spec[:5]) and r['prime']==prime and r['seed']==seed
            b,=r['blocks']
            assert b['Y']==spec[-1] and tuple(b[k] for k in ('dim','rank_in','rank_out','rho'))==expected
            assert b['rank_in']+b['rank_out']+b['rho']==b['dim']
            print(json.dumps(dict(prime=prime,seed=seed,sector=spec[:5],Y=b['Y'],dim=b['dim'],rank_I_lower=b['rank_in'],rank_joint_lower=b['rank_out']+b['rho'])))
