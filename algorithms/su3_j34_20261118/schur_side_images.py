"""Retain all side images when eliminating first-equation targets 1 and 5.

Rows refer to the original ordered source recipes. All arrays contain unsigned
little-endian residues. Tiling restricts the same evaluation points and hence
changes neither the field nor the matrix. A read-back audit is performed before
the final per-target receipt is written.
"""
import argparse
from copy import copy
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import schur_coordinates as sc

m = sc.m
maps = m.maps


def cut(frame, start, stop):
    """View a contiguous interval of points of an existing frame, without native caches."""
    f = copy(frame)
    for name in tuple(f.__dict__):
        if name.startswith('_native_su_'):
            delattr(f, name)
    f.P = stop - start
    f.ncols = f.cpp * f.P
    f.L = {l: (k, v[start:stop]) for l, (k, v) in frame.L.items()}
    return f


def pivot_images(a, t, words, p, width):
    base = a.sides / f'side_{t}'
    meta = json.loads((base/'result.json').read_bytes())
    manifest = json.loads((base/'manifest.json').read_bytes())
    active = [j for j, (i, _) in enumerate(words) if t in p['sources'][m.SOURCES[i]]['targets']]
    assert manifest['pivot_rows'] == active and meta['active_pivot_rows'] == len(active)
    assert meta['manifest_sha256'] == m.fsha(base/'manifest.json')
    assert meta['minor_sha256'] == manifest['minor_sha256'] == m.INPUT and manifest['plan_sha256'] == m.PLAN
    assert meta['width'] == manifest['width'] == width and meta['prime'] == manifest['prime'] == m.P
    assert meta['target'] == manifest['target'] == t and meta['entries'] == len(active)*width
    data = np.empty((len(active), width), dtype='<u4')
    h = hashlib.sha256()
    for j in range(0, len(active), manifest['batch']):
        block = np.load(base/f'{j}.npy', allow_pickle=False)
        assert block.shape == (min(manifest['batch'], len(active)-j),width) and block.dtype == data.dtype
        h.update(block.tobytes()); data[j:j+len(block)] = block
    assert h.hexdigest() == meta['entries_le_u4_sha256'] and np.count_nonzero(data) == meta['nonzero']
    return active, data, h.hexdigest()


def run(a):
    sys.path.insert(0, str(a.runtime.resolve()))
    from exact_blas_rows import matmul_mod
    from native_su_profile_derivation_v1 import evaluate
    p, words, _ = m.read_inputs(a.plan, a.minor)
    t = a.target
    maps.setup(a.runtime.resolve(), p, [t], list(m.SOURCES), m.P)
    f = maps.frames[0]; active, pv, pivot_hash = pivot_images(a, t, words, p, f.ncols)
    out = a.out / f'target_{t}'; out.mkdir(exist_ok=True, parents=True)
    # Sample one actual production point; no random frame reconstruction or reseeding.
    test_point = f.P//2; one = cut(f, test_point, test_point+1)
    sample = slice(test_point*f.cpp, (test_point+1)*f.cpp)
    tests = sorted({0, len(active)//2, len(active)-1})
    for j in tests:
        np.testing.assert_array_equal(sc.grouped_action(p, one, t, words[active[j]], maps.charge), pv[j,sample])
    records = []
    for s in m.SOURCES:
        idx, rows, nf, _, basis_hash = sc.source_basis(a, p, words, s)
        free = rows[nf:]; cert = json.loads((a.coords/f'coordinates_{s}.json').read_bytes())
        cf = a.coords/f'coefficients_{s}.npy'
        assert m.fsha(cf) == cert['coefficient_file_sha256'] and cert['source_basis_file_sha256'] == basis_hash
        c = np.load(cf, allow_pickle=False)
        assert c.shape == (len(free),len(words)) and c.dtype == np.dtype('<u4')
        assert cert['source_minor_determinant'] and cert['nonzero_residuals'] == 0 and cert['source'] == s and cert['pivot_minor_sha256'] == m.INPUT and cert['prime'] == m.P
        co = np.ascontiguousarray(c[:,active]); del c
        use = t in p['sources'][s]['targets']; allzero = not len(free) or (not np.any(co) and not use)
        rec = dict(source=s,free=len(free),width=f.ncols,source_basis_sha256=basis_hash,
                   coefficients_sha256=cert['coefficient_file_sha256'],free_reaches_target=use,
                   active_pivot_coefficients=int(np.count_nonzero(co)),structural_zero=allzero)
        if allzero:
            rec['nonzero_entries']=0;records.append(rec);continue
        path = out/f'source_{s}.npy'; temp = path.with_suffix('.part.npy')
        dst = np.lib.format.open_memmap(temp, mode='w+',dtype='<u4',shape=(len(free),f.ncols))
        points = [(k,min(k+a.tile_points,f.P)) for k in range(0,f.P,a.tile_points)]
        frames = [cut(f,k,l) for k,l in points] if use else []
        source = maps.pattern(p['sources'][s]); target = maps.pattern(p['targets'][t]); nnz = 0
        for off in range(0,len(free),a.rows):
            group = free[off:off+a.rows]; part = co[off:off+len(group)]
            rhs = (np.zeros((len(group),f.ncols),dtype=np.int64) if not part.any() else matmul_mod(part,pv,m.P))
            if use:
                expr = [maps.project(maps.algebra.parity_project({w:1},1),maps.charge,m.P) for w in group]
                for frame,(k,l) in zip(frames,points):
                    rhs[:,k*f.cpp:l*f.cpp] -= evaluate(expr,frame,source,target)
            v = (-rhs) % m.P; nnz += np.count_nonzero(v);dst[off:off+len(group)]=v
        dst.flush();del dst
        # Read the serialized rows; apply the ordinary integer dot product to all pivot
        # rows at this point, independently of BLAS. Native fresh-action checks expand
        # the sampled free words with the graded word differential, before evaluation.
        saved = np.load(temp,mmap_mode='r',allow_pickle=False); checks=[]
        for j in sorted({0,len(free)//2,len(free)-1}):
            v = (sc.grouped_action(p,one,t,(idx,free[j]),maps.charge) if use else np.zeros(f.cpp,dtype=np.int64))
            corr = [sum(int(coeff)*int(x) for coeff,x in zip(co[j],column)) % m.P for column in pv[:,sample].T]
            np.testing.assert_array_equal(saved[j,sample],(v-corr)%m.P);checks.append(j)
        assert np.count_nonzero(saved)==nnz; del saved; temp.replace(path)
        rec.update(path=path.name,sha256=m.fsha(path),nonzero_entries=int(nnz),direct_rows=checks,direct_point=test_point,
                   direct_coordinates=len(checks)*f.cpp);records.append(rec)
    result = dict(target=t,prime=m.P,minor_sha256=m.INPUT,plan_sha256=m.PLAN,pivot_entries_sha256=pivot_hash,
                  pivot_direct_rows=[active[j] for j in tests],point=test_point,width=f.ncols,rows=records,
                  implicit_or_stored_entries=sum(r['free']*f.ncols for r in records),sign='free image minus full pivot combination')
    m.save_json(out/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--runtime',type=Path,default=maps.RUNTIME)
    ap.add_argument('--plan',type=Path,default=maps.HERE/'classical_order1.json')
    ap.add_argument('--minor',type=Path,default=maps.HERE/'other_maps/task_5/minor.json.gz')
    for key in ('coords','sides','bases','out'):ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--target',type=int,choices=m.OTHERS,required=True);ap.add_argument('--tile-points',type=int,default=96);ap.add_argument('--rows',type=int,default=12)
    args=ap.parse_args();assert 0<args.tile_points<=1024 and 0<args.rows<=48;run(args)
