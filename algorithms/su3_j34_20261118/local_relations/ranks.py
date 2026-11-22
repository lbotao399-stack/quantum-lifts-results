"""Compute exact characters and joint rational lower ranks for the R20 repairs."""
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ALG = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ALG/'dual_singlets_20261118'))
from trace_descent import dimensions

SLOTS = (((3,2,2),1,6,(392,369),(312,426,23)),
         ((2,1,1),2,5,(225,205),(146,270,14)),
         ((1,0,0),3,4,(18,12),(7,21,2)),
         ((2,1,0),1,4,(8,4),(4,4,4)))


def verify():
    native = ALG/'dual_singlets_20261118/taylor_hooks/native'
    env = dict(os.environ,PAR_THREADS='2',A_CHECK='1',A_VERB='0',A_SCALAR='1')
    env.pop('A_NORHO',None)
    result = []
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'engine'
        subprocess.run(['c++','-O3','-std=c++17','-pthread',str(native/'engine.cpp'),'-o',str(exe)],check=True)
        for theta,t,length,expected_dims,expected_ranks in SLOTS:
            dim = tuple(sum(dimensions((t,t,*theta),length,3,c).values()) for c in ('SU+','SU-'))
            assert dim == expected_dims
            for prime,offset in ((1000003,0),(1000033,31)):
                seed = 28431 + (t-1)*93 + offset
                args = ['chain','SU3',*theta,t,t,length,length,seed,prime]
                row = json.loads(subprocess.check_output([str(exe),*map(str,args)],env=env))
                assert row['sector'] == [*theta,t,t] and row['prime'] == prime and row['seed'] == seed
                block, = row['blocks']
                ranks = tuple(block[k] for k in ('rank_in','rank_out','rho'))
                assert ranks == expected_ranks and sum(ranks) == sum(dim) == block['dim']
                assert block['Y'] == length and block['F'] == 1
            result.append(dict(theta=theta,dots=t,length=length,colour_dimensions=dim,
                               incoming=ranks[0],outgoing=ranks[1],Cartan_on_kernel=ranks[2]))
    return result


if __name__ == '__main__':
    if not __debug__:
        raise SystemExit('enable assertions')
    print(json.dumps(verify(),indent=2))
