"""Find pivot coefficients for each free trace direction, including independent checks on other targets."""
import argparse,gzip,hashlib,json,math,sys
from pathlib import Path
import numpy as np
from flint import nmod_mat
import schur_first as m
maps=m.maps

def source_basis(a,p,words,s):
    path=a.bases/f'basis_{s}.json.gz';record=json.loads((a.bases/f'basis_{s}.json').read_bytes());assert record['complete'] and record['sha256']==m.fsha(path) and record['minor_sha256']==m.INPUT
    b=json.loads(gzip.decompress(path.read_bytes()));rows=[maps.decode(w) for w in b['words']];i=m.SOURCES.index(s);fixed=[w for k,w in words if k==i];d=p['sources'][s]['dimension'];assert b['source']==s and len(rows)==len(set(rows))==record['rank']==d==b['dimension'] and rows[:len(fixed)]==fixed and len(fixed)==b['fixed'];cols=b['columns'];assert len(cols)==len(set(cols))==d
    return i,rows,len(fixed),cols,record['sha256']

def grouped_action(p,f,s,raw,charge):
    # Separate implementation: expand the exact graded word differential before evaluating emitted traces.
    from native_su_values_v1 import evaluate
    from sp_full_profiles import profile
    i,w=raw;src=m.SOURCES[i];expr=maps.project(maps.algebra.parity_project({w:1},1),charge,m.P);out=maps.algebra.differential(expr);out={w:int(c)%m.P for w,c in out.items() if c%m.P and profile(tuple(l for t in w for l in t))==maps.pattern(p['targets'][s])}
    if s not in p['sources'][src]['targets']:assert not out
    return evaluate([out],f)[0]

def direct(a,p,words,cols,s,free,co,out):
    from fortuity.full.evalg import Frame
    from g2_full_classical_screen_v1 import contents,sector
    from exact_blas_rows import matmul_mod
    alpha=sorted({l for n in (p['source_length'],p['target_length']) for c in contents(tuple(p['mu']),n) for l in c});charge=sector(tuple(p['mu']),p['source_length']);samples=sorted(set((0,len(free)//2,len(free)-1)));sel=[free[i] for i in samples];results=[];arrays={};assert all(i>=0 for i in samples)
    for t,row in enumerate(p['targets']):
        d=row['fermions'];f=Frame(3,1,d+2,[d],alpha,seed=2026111941+7919*t,p=m.P);width=f.ncols
        # Vectors on all 17 profiles. Structural zeros are left explicit. For all other
        # pivot directions compute d in the pinned word algebra, with no block splitting.
        pv=np.zeros((len(words),width),dtype=np.int64)
        if any(t in p['sources'][m.SOURCES[k]]['targets'] for k,w in words):
            for j,r in enumerate(words):
                if t in p['sources'][m.SOURCES[r[0]]]['targets']:pv[j]=grouped_action(p,f,t,r,charge)
        fv=np.vstack([grouped_action(p,f,t,r,charge) if t in p['sources'][m.SOURCES[r[0]]]['targets'] else np.zeros(width,dtype=np.int64) for r in sel]);delta=(fv-matmul_mod(co[samples].astype(np.int64),pv,m.P))%m.P
        if t in m.TARGETS:assert not delta.any()
        if t not in m.TARGETS+m.OTHERS:assert not delta.any()
        arrays[f't{t}']=m.pack(delta);results.append(dict(target=t,width=width,nonzero=int(np.count_nonzero(delta)),sha256=m.ahash(delta)))
    np.savez_compressed(out/f'fresh_{s}.npz',**arrays);return dict(profile=s,free_indices=samples,points_per_target=1,method='direct word differential before target specialization',targets=results,archive_sha256=m.fsha(out/f'fresh_{s}.npz'))

def run(a):
    sys.path.insert(0,str(a.runtime.resolve()));from exact_blas_rows import matmul_mod
    from native_su_parallel_v1 import evaluate
    p,words,cols=m.read_inputs(a.plan,a.minor);s=a.profile;i,rows,k,src_cols,basis_hash=source_basis(a,p,words,s);free=[(i,w) for w in rows[k:]];out=a.out;out.mkdir(exist_ok=True,parents=True)
    maps.setup(a.runtime.resolve(),p,list(m.TARGETS),list(m.SOURCES),m.P)
    factor=json.loads((a.factor/'factor.json').read_bytes());inverse=a.factor/'target_1_5_inverse.npy';assert factor['minor_sha256']==m.INPUT and m.fsha(inverse)==factor['inverse_file_sha256'];inv=np.load(inverse,allow_pickle=False);assert inv.shape==(1112,1112)
    # Replay the full-rank source minor independently of the basis selection's incremental reduction.
    f=m.basis_frame(p,s);M=nmod_mat(len(rows),len(rows),m.P)
    for j in range(0,len(rows),24):
        b=rows[j:j+24];x=[maps.project(maps.algebra.parity_project({w:1},1),maps.charge,m.P) for w in b];v=evaluate(x,f,a.workers)[:,src_cols]
        for n,one in enumerate(v):
            for l,val in enumerate(one.tolist()):M[j+n,l]=val
    det=int(M.det());assert det;del M
    F=np.vstack([maps.values(g) for g in m.batches(free,12)]) if free else np.zeros((0,1131),dtype=np.int64);B=np.vstack([maps.values(g) for g in m.batches(words,12)]);assert m.ahash(B)==factor['full_projected_entries_sha256']
    co=matmul_mod(F[:,cols],inv,m.P) if free else np.zeros((0,1112),dtype=np.int64);check=matmul_mod(co,B,m.P) if free else np.zeros_like(F);assert np.array_equal(check,F);co=m.pack(co);path=out/f'coefficients_{s}.npy';np.save(path,co,allow_pickle=False)
    r=dict(source=s,source_recipes=len(rows),source_basis_file_sha256=basis_hash,source_minor_order=len(rows),source_minor_determinant=det,pivots=k,free=len(free),target_coords=1131,coordinates_checked=len(free)*1131,nonzero_residuals=0,nonzero_coefficient_entries=int(np.count_nonzero(co)),coefficient_file_sha256=m.fsha(path),coefficient_entries_sha256=m.ahash(co),pivot_minor_sha256=m.INPUT,prime=m.P,coefficient_sign='subtract stored combination of full pivot recipes from each free recipe');m.save_json(out/f'coordinates_{s}.json',r);print(json.dumps(r),flush=True)
    if a.fresh and free:m.save_json(out/f'fresh_{s}.json',direct(a,p,words,cols,s,free,co,out))

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,default=maps.HERE/'classical_order1.json');ap.add_argument('--minor',type=Path,default=maps.HERE/'other_maps/task_5/minor.json.gz');ap.add_argument('--runtime',type=Path,default=maps.RUNTIME);ap.add_argument('--bases',type=Path,required=True);ap.add_argument('--factor',type=Path,required=True);ap.add_argument('--profile',type=int,choices=m.SOURCES,required=True);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--fresh',action='store_true');run(ap.parse_args())
