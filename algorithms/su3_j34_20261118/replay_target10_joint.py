"""Independent exact replay of the saved selected-row solve; no generating code is imported."""
import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from flint import nmod_mat


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(values):
    return hashlib.sha256(np.asarray(values, dtype='<u4').tobytes()).hexdigest()


def run(folder, root):
    r = json.loads((folder/'result.json').read_bytes())
    p = r['prime']; assert p == 1000003
    files = {'recipes.json.gz':'recipe_sha256', 'matrix.npy':'matrix_sha256', 'solution.npz':'solution_sha256'}
    hashes = {name:sha(folder/name) for name in files}
    assert all(hashes[name] == r[key] for name,key in files.items())
    assert sha(root/'classical_order1.json') == r['plan_sha256']
    plan = json.loads((root/'classical_order1.json').read_bytes())['incoming']
    words = json.loads(gzip.decompress((folder/'recipes.json.gz').read_bytes()))
    supports = {s: set(plan['sources'][s]['targets']) for s,_ in words}
    assert set(supports) == set(r['sources'])
    active = sorted(set.union(*supports.values())); assert active == r['reachable']
    assert {str(k):v for k,v in Counter(s for s,_ in words).items()} == r['input_recipes_by_source']
    M = np.load(folder/'matrix.npy', allow_pickle=False)
    with np.load(folder/'solution.npz', allow_pickle=False) as a:
        c, b, residual, cols = [a[k].tolist() for k in ('coefficients','rhs','residual','columns')]
    assert M.shape == (len(words),len(b)) == (704,2244) and len(c)==704
    # Arbitrary-precision Python integers; each dot product runs in the opposite orientation from the producer.
    got = [(y - sum(x*int(z) for x,z in zip(c,col))) % p for y,col in zip(b,M.T)]
    assert got == residual
    stats, train, hold, off = {}, [], [], 0
    for i,t in enumerate(plan['targets']):
        assert t['index']==i
        from math import comb
        width = comb(t['fermions']+2,t['fermions']); cut,end = off+8*width,off+12*width
        a, h = got[off:cut], got[cut:end]; old = r['all_target_stats'][str(i)]
        assert digest(a+h)==old['residual_sha256'] and sum(map(bool,a))==old['training_nonzero'] and sum(map(bool,h))==old['withheld_nonzero']
        if i in (6,10,13,14,15,16):
            f = (root/('joint/edges' if i in (6,10) else 'first_profile16/coupled')/f'fresh_{i}.npy')
            assert sha(f)==r['baseline_hashes'][str(i)]
            samples = np.load(f,allow_pickle=False); rhs = samples if i in (6,10) else samples[2]
            assert rhs.tolist()==b[off:end]
        else: assert not any(b[off:end])
        ix=[k for k,(s,_) in enumerate(words) if i not in supports[s]]
        assert not np.any(M[ix,off:end]); stats[str(i)]={'training_nonzero':sum(map(bool,a)), 'withheld_nonzero':sum(map(bool,h)), 'coordinates':end-off}
        if i in active: train+=list(range(off,cut)); hold+=list(range(cut,end))
        off=end
    assert off==len(got) and set(train)==set(cols) and len(cols)==len(set(cols))==704
    T = M[:,cols]; det = int(nmod_mat(T.tolist(),p).det()); assert det==r['minor_determinant'] and det and digest(T)==r['minor_sha256']
    # One augmented minor with the original RHS as the last row witnesses incompatibility for the entire selected span.
    j=next(k for k in hold if got[k]); augmented=np.vstack((M[:,cols+[j]],np.asarray([b[k] for k in cols+[j]],dtype='<u4')))
    witness=int(nmod_mat(augmented.tolist(),p).det()); assert witness==det*got[j]%p and witness
    return dict(prime=p,files=hashes,matrix_shape=list(M.shape),all_residual_coordinates=len(got),all_nonzero=sum(map(bool,got)),training_coordinates=len(train),training_nonzero=sum(bool(got[k]) for k in train),withheld_coordinates=len(hold),withheld_nonzero=sum(bool(got[k]) for k in hold),selected_minor_order=len(cols),selected_minor_determinant=det,augmented_minor_order=len(cols)+1,extra_column_zero_based=j,extra_residual=got[j],augmented_minor_determinant=witness,per_target=stats)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('folder',type=Path);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--output',type=Path);args=ap.parse_args()
    assert __debug__; res=run(args.folder,args.root);text=json.dumps(res,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')
