"""Replay source/target assembly, both ranks, pivot determinant and the solution without BLAS.

The input is span.py's certificate and independent NPY blocks. No solver functions
or pivot selection from span.py are invoked. Both ranks use FLINT on the entire matrix.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from flint import nmod_mat


def digest(x):
    return hashlib.sha256(x).hexdigest()


def load(f, pin, dims, p):
    assert digest(f.read_bytes()) == pin, f
    x = np.load(f, allow_pickle=False)
    assert x.shape == dims and x.dtype == np.dtype('<u4') and x.max(initial=0) < p, f
    return x


def verify(op, cert, blocks, rhs):
    p = cert['prime']; ts = cert['target_profiles']; ss = cert['sources']
    plan = Path(__file__).with_name(op + '_plan.json')
    assert digest(plan.read_bytes()) == cert['plan_sha256']; plan = json.loads(plan.read_bytes())
    js = [t['index'] for t in ts]; assert len(set(js)) == len(js)
    assert [s['index'] for s in ss] == [s['index'] for s in plan['sources'] if set(s['targets']) & set(js)]
    v = np.zeros((cert['source_rows'], cert['coordinates']), dtype='<u4'); bs = []; n = 0
    for source in ss:
        i, d = source['index'], source['dimension']; assert d == plan['sources'][i]['dimension']; col = 0
        assert set(map(int,source['blocks'])) == set(js) & set(plan['sources'][i]['targets'])
        for target in ts:
            j, w = target['index'], target['columns']; assert (i in plan['targets'][j]['sources']) == (str(j) in source['blocks'])
            if str(j) in source['blocks']:
                f = blocks / op / f'sources{i:02d}' / f't{j:02d}.npy'
                v[n:n+d, col:col+w] = load(f, source['blocks'][str(j)], (d,w), p)
            col += w
        assert col == v.shape[1]; n += d
    assert n == v.shape[0]
    for target in ts:
        j, w = target['index'], target['columns']; assert target['dimension'] == plan['targets'][j]['dimension']
        bs.append(load(rhs / f'{op}_t{j:02d}.npy',target['values_sha256'],(1,w),p))
    b = np.hstack(bs); assert digest(v.tobytes()) == cert['matrix_row_major_sha256'] and digest(b.tobytes()) == cert['rhs_row_major_sha256']
    def as_mat(a):return nmod_mat(a.tolist(),p)
    rank = as_mat(v).rank(); aug = as_mat(np.vstack((v,b))).rank()
    rows, cols = cert['minor_rows'], cert['minor_columns']; assert rank == cert['rank'] == len(rows) == len(cols) and aug == cert['augmented_rank']
    det = int(as_mat(v[np.ix_(rows, cols)]).det()); assert det == cert['minor_determinant'] != 0
    if rank == aug:
        co = cert['solution_on_minor_rows']; assert len(co) == len(rows) and all(0 <= z < p for z in co)
        # Each unreduced coordinate is bounded by rank*(p-1)^2, including every partial sum.
        assert rank*(p-1)**2 < 2**64
        product = (np.array(co,dtype=np.uint64) @ v[rows].astype(np.uint64)) % p
        assert np.array_equal(product,b[0]); residuals = int(np.count_nonzero(product != b[0])); pairing = None
    else:
        f = cert['separator']; assert len(f)==v.shape[1] and all(0 <= z < p for z in f) and len(f)*(p-1)**2 < 2**64
        f = np.array(f,dtype=np.uint64); residuals = sum(int(np.count_nonzero((x.astype(np.uint64) @ f)%p)) for x in np.array_split(v,16)); pairing = int((b.astype(np.uint64) @ f)[0])%p
        assert residuals == 0 and pairing == cert['pairing'] != 0
    return dict(prime=p,profiles=js,source_rows=n,coordinates=v.shape[1],rank=rank,augmented_rank=aug,minor_determinant=det,nonzero_residuals=residuals,separator_pairing=pairing)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ('certificate','blocks','rhs','output'):ap.add_argument(k,type=Path)
    x=ap.parse_args();assert __debug__;cert=json.loads(x.certificate.read_bytes());out={}
    for op,v in cert['results'].items():
        out[op]=verify(op,v,x.blocks,x.rhs);print(op,out[op],flush=True)
    result=dict(certificate_sha256=digest(x.certificate.read_bytes()),algorithm='Direct FLINT ranks of complete rectangular matrices; bounded uint64 integer residuals.',results=out)
    tmp=x.output.with_suffix('.tmp');tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(x.output)
