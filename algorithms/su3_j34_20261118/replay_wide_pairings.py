"""Replay the saved profile-10 separator with bounded integer operations, without the solver."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

PRIME = 1000003
RECORD = 'd39a6316ed21d22281a3b62b12526f93ee338ba4363762fe402aba9145aa0c56'
SEPARATOR = 'db15ae166ff3d7f4170f84998aa74a6bd3bc7a6c23b67b7a868eb6ceb4699807'


def file_hash(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def run(arrays, certificate):
    files = {n: file_hash(p/n) for p, n in ((arrays,'matrix.npy'), (arrays,'rhs.npy'),
             (certificate,'wide.json'), (certificate,'separator.npy'))}
    assert files['wide.json'] == RECORD and files['separator.npy'] == SEPARATOR
    record = json.loads((certificate/'wide.json').read_bytes())
    matrix, rhs, dual = (np.load(p/n, mmap_mode='r', allow_pickle=False)
        for p,n in ((arrays,'matrix.npy'),(arrays,'rhs.npy'),(certificate,'separator.npy')))
    assert (matrix.shape,rhs.shape,dual.shape) == ((6848,7000),(7000,),(7000,))
    assert all(x.dtype.str == '<u4' and x.flags.c_contiguous for x in (matrix,rhs,dual))
    assert record['prime'] == PRIME and max(rhs.max(),dual.max()) < PRIME
    assert len(rhs)*(PRIME-1)**2 < np.iinfo(np.int64).max
    raw = hashlib.sha256(); evaluations = hashlib.sha256(); rows=0; nonzero=0
    dual64 = dual.astype(np.int64)
    for start in range(0,len(matrix),37):
        block = matrix[start:start+37]; raw.update(memoryview(block).cast('B'))
        assert block.max() < PRIME
        dots = (block.astype(np.int64) @ dual64) % PRIME
        nonzero += int(np.count_nonzero(dots)); rows += len(dots)
        evaluations.update(dots.astype('<u4').tobytes())
    right_hash = hashlib.sha256(rhs.tobytes()).hexdigest()
    assert (raw.hexdigest(),right_hash) == (record['matrix_sha256'],record['rhs_sha256'])
    pair = sum(map(lambda xy: int(xy[0])*int(xy[1]),zip(rhs,dual))) % PRIME
    assert rows == 6848 and nonzero == 0 and pair == record['pairing'] == 923594
    n = len(record['pivot_rows']); assert n == len(set(record['pivot_rows'])) == len(set(record['pivot_columns'])) == record['rank'] == 6635
    assert math.prod(record['panel_determinants']) % PRIME == record['minor_determinant'] == 650555
    assert pair * record['minor_determinant'] % PRIME == record['augmented_minor_determinant'] == 892132
    return dict(prime=PRIME,rows_replayed=rows,width=len(rhs),nonzero_row_pairings=nonzero,rhs_pairing=pair,
                file_sha256=files,data_sha256={'matrix':raw.hexdigest(),'rhs':right_hash,'row_pairings':evaluations.hexdigest()},
                feed='original archived arrays; 37-row integer blocks; Python integer RHS',
                minor_check='recorded panel determinant product and Schur identity only; no independent factorization here')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument('arrays',type=Path); ap.add_argument('--certificate',type=Path)
    a=ap.parse_args(); print(json.dumps(run(a.arrays,a.certificate or a.arrays),indent=2))
