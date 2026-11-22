"""Solve the first correction on all evaluation coordinates and test the second residual."""
from argparse import ArgumentParser
from pathlib import Path
import hashlib
import json
import struct
import numpy as np
from flint import nmod_mat


def read(path):
    data = np.fromfile(path, dtype='<u4')
    magic, p, h, w = map(int, data[:4])
    assert magic == 0x47325130 and len(data) == 4+h*w
    return p, data[4:].reshape(h,w)


def write(path, prime, array):
    contiguous = np.ascontiguousarray(array, dtype='<u4')
    path.write_bytes(struct.pack('<IIII',0x47325130,prime,*contiguous.shape)+contiguous.tobytes())


def matrix(array, p):
    return nmod_mat(array.tolist(),p)


def run(directory):
    def get(stem):
        p, raw = read(directory/(stem+'.u32'))
        assert p == prime
        return raw
    prime, boundary = read(directory/'in_1_9.u32')
    raw = get('acted_cycle_0_9')
    rhs = np.where(raw==0,0,prime-raw)
    pivots = list(map(int,(directory/'pivots_1_8.u32').read_text().split()))
    columns = boundary[:,pivots]
    echelon, rank = matrix(columns.T,prime).rref()
    assert rank == len(pivots)
    selected_rows = np.asarray(echelon.tolist(),dtype=np.uint32).astype(bool).argmax(axis=1)
    square = matrix(columns[selected_rows],prime)
    det = int(square.det())
    assert det
    coeff = square.solve(matrix(rhs[selected_rows],prime))
    assert matrix(columns,prime)*coeff == matrix(rhs,prime)
    correction = np.zeros((raw.shape[1],boundary.shape[1]), dtype=np.uint32)
    correction[:,pivots]=np.asarray(coeff.tolist(),dtype=np.uint32).T
    write(directory/'primitive_1_8.u32',prime,correction)
    first = get('paired_cycle_1_8')
    result = dict(prime=prime,first_image_rank=rank,minor=det,
                  first_source_vectors=raw.shape[1],first_residual_coordinates=rhs.size,
                  bare_first_nonzero=int(np.count_nonzero(raw)),
                  homogeneous_corrections=first.shape[1],
                  homogeneous_first_pairing=first.tolist(),
                  correction_terms=[int(np.count_nonzero(row)) for row in correction])
    if (directory/'acted_primitive_1_8.u32').exists():
        raw_second=get('acted_primitive_1_8')
        state = matrix(get('eval_2_8'),prime).solve(matrix(raw_second,prime))
        closed=matrix(get('in_2_9'),prime)*state
        assert all(not x for row in closed.tolist() for x in row)
        projected=matrix(get('dual_2_8'),prime)*matrix(raw_second,prime)
        result.update(second_residual_coordinates=raw_second.size,
                      raw_second_nonzero=int(np.count_nonzero(raw_second)),
                      second_residual_classical_closure=[[int(x) for x in row] for row in closed.tolist()],
                      second_pairing=[[int(x) for x in row] for row in projected.tolist()])
        if not np.any(first) and not any(x for row in projected.tolist() for x in row):
            # Exhibit an actual second primitive on *all* target coordinates.
            before=get('in_2_8')
            independent=list(map(int,(directory/'pivots_2_7.u32').read_text().split())) if (directory/'pivots_2_7.u32').exists() else None
            if independent is None:
                echelon, r=matrix(before,prime).rref()
                independent=np.asarray(echelon.tolist(),dtype=np.uint32)[:r].astype(bool).argmax(axis=1)
            selected=before[:,independent]
            te, r=matrix(selected.T,prime).rref()
            pts=np.asarray(te.tolist(),dtype=np.uint32)[:r].astype(bool).argmax(axis=1)
            new=matrix(selected[pts],prime).solve(-matrix(raw_second[pts],prime))
            assert matrix(selected,prime)*new==-matrix(raw_second,prime)
            primitive=np.zeros((raw_second.shape[1],before.shape[1]),dtype=np.uint32)
            primitive[:,independent]=np.asarray(new.tolist(),dtype=np.uint32).T
            write(directory/'primitive_2_7.u32',prime,primitive)
            result['second_equation_coordinates']=raw_second.size
            result['second_primitive_terms']=[int(np.count_nonzero(row)) for row in primitive]
    files = sorted(directory.glob('*.u32'))
    result['sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    (directory/'obstruction.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ('sha256','homogeneous_first_pairing','second_residual_classical_closure')},flush=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    if not __debug__: raise RuntimeError('enable assertions')
    run(p.parse_args().directory)
