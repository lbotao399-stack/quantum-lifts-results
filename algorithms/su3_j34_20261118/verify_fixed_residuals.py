"""Independently verify the stored first-correction residuals, not membership in the full image."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
P = 1000003

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check(root=HERE):
    folder = root/'first_profile16'
    sol, minor = folder/'solution.json.gz', root/'other_maps/task_3/minor.json.gz'
    solution = json.loads(gzip.decompress(sol.read_bytes()))
    cert = json.loads(gzip.decompress(minor.read_bytes()))
    projection = json.loads((folder/'profile_result.json').read_bytes())
    assert cert['map']['sources'] == [19,20,21,22] and projection['initial_minor_determinant'] == cert['determinant'] == 258962
    assert len(solution['base_coefficients']) == len(cert['rows']) == 7211 and all(0 < i < P for i in solution['base_coefficients'])
    assert not solution['extra_coefficients'] and not solution['extra_rows'] and solution['prime'] == projection['prime'] == P
    assert digest(sol) == projection['solution_sha256'] and digest(minor) == solution['base_minor_sha256']
    assert projection['projected_equation_solved'] and projection['replayed_coordinates'] == 7272 and projection['residual_nonzero'] == 0
    plan = root/'classical_order1.json'; facts = json.loads(plan.read_bytes())['incoming']; report = {}
    for mode,t in [('fresh',i) for i in (13,14,15,16)]+[('production',i) for i in (13,14,15)]:
        stem = f'{mode}_{t}'; file = folder/'coupled'/(stem+'.npy'); rec = json.loads(file.with_suffix('.json').read_bytes())
        assert (rec['mode'],rec['target'],rec['prime']) == (mode,t,P)
        assert rec['source_minor_sha256'] == digest(minor) and rec['solution_sha256'] == digest(sol) and rec['plan_sha256'] == digest(plan)
        assert digest(file) == rec['vectors_sha256'] and rec['target_dimension'] == facts['targets'][t]['dimension']
        a = np.load(file, allow_pickle=False); assert a.dtype.str == '<u4' and a.shape == (3, rec['width']) and np.all(a < P)
        x = a.astype(np.int64); assert np.array_equal((x[0]-x[1]) % P, x[2])
        for row,name in zip(a,('rhs','lhs','residual')):
            assert hashlib.sha256(row.tobytes()).hexdigest() == rec[name+'_sha256'] and np.count_nonzero(row) == rec[name+'_nonzero']
        nz = np.flatnonzero(x[2]); j = int(nz[0]) if len(nz) else None
        assert rec['first_nonzero'] == (None if j is None else dict(index=j,rhs=int(x[0,j]),lhs=int(x[1,j]),residual=int(x[2,j])))
        report[stem] = dict(coordinates=int(a.shape[1]),nonzero_residuals=len(nz))
    total = sum(r['coordinates'] for r in report.values()); bad = sum(r['nonzero_residuals'] for r in report.values())
    assert (total,bad)==(84894,84462)
    return dict(prime=P,entrywise_replay=report,coordinates=total,nonzero_residuals=bad,original_profile16_production_receipt={'coordinates':7272,'nonzero_residuals':0},scope='Independent signed subtraction and hashes of downloaded arrays. Original profile-16 production result read from its separate solve receipt; the full joint image is not tested.')

if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path);a=p.parse_args();text=json.dumps(check(),indent=2)+'\n'
    if a.out:a.out.write_text(text)
    print(text,end='')
