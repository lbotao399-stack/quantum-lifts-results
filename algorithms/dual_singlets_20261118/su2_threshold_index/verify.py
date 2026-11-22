"""Check every opposite-statistics rung in the SU(2) threshold index coefficient."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

from abelian import check


HERE = Path(__file__).resolve().parent
NATIVE = HERE.parent / 'taylor_hooks' / 'native'
FIELDS = ('Y', 'F', 'dim', 'rank_in', 'rank_out', 'h', 'rho')


def chain(engine, degrees, first, last, seed, prime):
    args = [str(engine), 'chain', 'SU2', *map(str, (*degrees, first, last, seed, prime))]
    env = dict(os.environ, PAR_THREADS='2', A_VERB='0', A_CHECK='1')
    env.pop('A_NORHO', None)
    raw = json.loads(subprocess.check_output(args, env=env))
    assert raw['gauge'] == 'SU2' and raw['sector'] == list(degrees)
    assert raw['seed'] == seed and raw['prime'] == prime
    blocks = [{key:b.get(key, 0) for key in FIELDS} for b in raw['blocks']]
    assert [b['Y'] for b in blocks] == list(range(first, last+1))
    for b in blocks:
        assert b['dim'] == b['rank_in'] + b['rank_out'] + b['h']
    return blocks


def verify(engine):
    columns = []
    for j, prime in enumerate((2147483647, 1000003)):
        primary_test = chain(engine, (4,4,3,0,0), 7, 7, 780541+97*j, prime)
        assert (primary_test[0]['dim'], primary_test[0]['rank_in'], primary_test[0]['rank_out']) == (11,6,5)
        opposite = []
        for t in range(5):
            degrees = (4-t,4-t,4-t,t,t)
            blocks = chain(engine, degrees, 1, sum(degrees), 783775+97*j, prime)
            if t == 0:
                core = blocks[6]
                assert (core['dim'], core['rank_in'], core['rank_out'], core['h']) == (11,4,6,1)
            for b in blocks:
                rung = b['Y'] + t - 7
                if rung % 2:
                    assert b['F'] == 0 and b['h'] == b['rho']
                    opposite.append(dict(t=t, rung=rung, **b))
        assert len(opposite) == 26
        columns.append(opposite)
    assert columns[0] == columns[1]
    survivors = [(b['t'],b['Y'],b['h']) for b in columns[0] if b['h']]
    assert survivors == [(2,4,3),(2,6,11),(4,4,1)]
    minors = [check(t, length) for t, length, _ in survivors]
    assert [(m['columns'],m['determinant']) for m in minors] == [(737,130460),(711,623010),(11,993091)]
    return dict(classical_primes=[2147483647,1000003], opposite_blocks=columns[0],
                nonzero_modular_cohomology=survivors, abelian_lowest_weight=minors)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine', type=Path, help='reuse a build of the pinned engine')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('enable assertions')
    for file, digest in json.loads((NATIVE/'sha256.json').read_bytes()).items():
        assert hashlib.sha256((NATIVE/file).read_bytes()).hexdigest() == digest
    with tempfile.TemporaryDirectory() as tmp:
        engine = args.engine.resolve() if args.engine else Path(tmp)/'engine'
        if args.engine is None:
            subprocess.run(['c++','-O3','-std=c++17','-pthread',str(NATIVE/'engine.cpp'),'-o',str(engine)], check=True)
        result = verify(engine)
    if args.out:
        args.out.write_text(json.dumps(result, indent=2)+'\n')
    print('26 opposite-parity blocks, 3 nonempty modular quotients, 0 compatible lowest weights')
