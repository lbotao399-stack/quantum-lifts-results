"""Test every derivative-free, nine-letter level-26 source of an outgoing d2."""
import argparse
import json
import os
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parent
PRIME = 2147483647
SOURCES = ((9,2,2),(8,3,2),(7,3,3),(7,4,2),(5,4,4),(6,5,2),(5,5,3),(6,4,3))


def block(engine, out, weight, length, seed, workers):
    key = '-'.join(map(str,weight)) + f'-Y{length}'
    path = out / (key + '.json')
    if not path.exists():
        args = [engine, *weight[:3], length, length, seed, workers, PRIME,
                '--probe-last', *weight[3:]]
        with (out / (key + '.log')).open('w') as log:
            raw = subprocess.check_output(list(map(str,args)), stderr=log, timeout=1500)
        record = json.loads(raw)
        assert record['weight'] == list(weight) and record['seed'] == seed
        path.write_text(json.dumps(record,indent=2)+'\n')
    record = json.loads(path.read_bytes())
    row, = record['blocks']
    assert record['weight'] == list(weight) and row['length'] == length and record['prime']==PRIME
    assert row['incoming']+row['outgoing']<=row['dimension']
    h = row['dimension']-row['incoming']-row['outgoing']
    print(weight, length, row['dimension'],row['incoming'],row['outgoing'],h,flush=True)
    return h


def search(engine, out, workers):
    out.mkdir(parents=True,exist_ok=True)
    assert set(SOURCES) == {(13-b-c,b,c) for b in range(2,10) for c in range(2,b+1)
                            if 13-b-c>=b}
    for i,w in enumerate(SOURCES):
        seed = 930371+601*i
        if block(engine,out,(*w,0,0),9,seed,workers):
            # A positive remainder is only an upper bound. Determine target
            # feasibility before constructing a rational source or quantum maps.
            block(engine,out,(*(v-2 for v in w),2,2),8,seed+149,workers)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--engine',type=Path,default=ROOT/'g2-native/g2_dotted_chain')
    p.add_argument('--out',type=Path,default=ROOT/'g2-level26-bmn')
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--lease',action='store_true')
    a=p.parse_args()
    if not __debug__ or a.workers<1: raise ValueError('positive workers and assertions required')
    bus=Path('~/.agent-bus/bin/agent-bus').expanduser()
    if a.lease:
        subprocess.run([bus,'lease',str(os.getpid()),'--agent','codex','--purpose',
                        'N=4 G2 level-26 tree and target rank tests','--hours','4'],check=True)
    try:
        search(a.engine.resolve(),a.out,a.workers)
    finally:
        if a.lease: subprocess.run([bus,'release',str(os.getpid())],check=True)
