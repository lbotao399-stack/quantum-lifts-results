"""Solve the level-34 first correction on profile 16, preserving the untested equations."""
import argparse, gzip, hashlib, json, math, sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import nmod_mat
import projected_q0 as maps

P=1000003
PLAN_HASH='1b7ac3dec77fcb50473d01176a8f309773e8411f0fd4e969f58299a026df6e8c'
MINOR_HASH='6db12024f5a86596a8c161628e0661097dc3b11b9cbee49019e9b9ce9542f199'
PROFILE_HASH='19186d93069a07ad04b3b29cf64664137d860d51e75e0bddefb159ef8f98c563'
RHS_HASH='9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14'

def sha(file):
    with file.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def pack_rhs(file,out):
    from sp_full_profiles import profile
    target=maps.pattern(json.loads((maps.HERE/'classical_order1.json').read_text())['incoming']['targets'][16])
    h=hashlib.sha256();n=terms=0
    with gzip.open(file,'rb') as f, out.open('wb') as disk,gzip.GzipFile(fileobj=disk,mode='wb',mtime=0,filename='') as dest:
        for line in f:
            h.update(line); c,w=json.loads(line); w=maps.decode(w);terms+=1
            if profile(tuple(l for tr in w for l in tr))==target:
                dest.write(line);n+=1
    assert h.hexdigest()==RHS_HASH and terms==517468 and n==87234
    print(json.dumps(dict(profile=16,rows=n,sha256=sha(out))),flush=True)

def in_field(c):
    x=Fraction(c)
    return (x.numerator*pow(x.denominator,-1,P))%P

def setup(runtime,plan,targets,sources,prime):
    maps.setup(runtime,plan,targets,sources,prime)

def rhs_values(args):
    from native_su_parallel_v1 import evaluate
    from sp_full_profiles import profile
    frame=maps.frames[0]; pattern=maps.pattern(maps.target_rows[0])
    result=np.zeros(frame.ncols,dtype=np.int64)
    expr={};n=0;checked=False
    def add(expr,checked):
        nonlocal result
        pairs=list(expr.items())
        if not checked and any(Fraction(c).denominator>1 for _,c in pairs):
            samples=[(w,c) for w,c in pairs if Fraction(c).denominator>1][:3]
            # Independently evaluate on the first full point, multiplying by a common
            # denominator before reducing. This includes every exterior coordinate there.
            import copy
            small=copy.copy(frame);small.P=1;small.ncols=frame.cpp
            small.L={l:(d,x[:1]) for l,(d,x) in frame.L.items()}
            base=small.eval_words([w for w,c in samples],dtype=np.int64)
            den=math.lcm(*(Fraction(c).denominator for _,c in samples))
            weights=np.array([int(Fraction(c)*den)%P for _,c in samples])
            ref=(weights@base%P)*pow(den,-1,P)%P
            reduced={w:in_field(c) for w,c in samples}
            assert np.array_equal(evaluate([reduced],frame,args.workers)[0,:frame.cpp],ref)
            checked=True
        result=(result+evaluate([{w:in_field(c) for w,c in pairs}],frame,args.workers)[0])%P
        return checked
    with gzip.open(args.rhs,'rb') as f:
        for line in f:
            c,w=json.loads(line);w=maps.decode(w)
            assert profile(tuple(l for t in w for l in t))==pattern and w not in expr
            expr[w]=c;n+=1
            if len(expr)==512:
                checked=add(expr,checked);expr={}
        if expr:checked=add(expr,checked)
    assert n==87234 and checked and len(result)==7272
    np.save(args.output/'rhs.npy',result.astype('<u4'),allow_pickle=False)
    print(json.dumps(dict(terms=n,width=len(result),nonzero=int(np.count_nonzero(result)),sha256=sha(args.output/'rhs.npy'))),flush=True)

def solve(args,cert):
    from exact_blas_rows import matmul_mod
    from su_modular_rowspace_v1 import RowSpace
    plan=json.loads(args.plan.read_text())['incoming'];sources=cert['map']['sources']; cols=cert['columns']; n=len(cols)
    rows=[(i,maps.decode(w)) for i,w in cert['rows']]
    init=(args.runtime.resolve(),plan,[16],sources,P)
    assert n==7211 and len(set(cols))==n and len(rows)==n and cert['determinant']==258962
    frame=maps.frames[0];width=frame.ncols;assert width==7272
    rhs=np.load(args.output/'rhs.npy',allow_pickle=False).astype(np.int64);assert rhs.shape==(width,)
    disk=args.output/'boundary.npy'
    with ProcessPoolExecutor(max_workers=args.workers,initializer=setup,initargs=init) as pool:
        if not disk.exists():
            dest=np.lib.format.open_memmap(args.output/'boundary.tmp',mode='w+',dtype='<u4',shape=(n,width))
            for i,packet in enumerate(pool.map(maps.values,(rows[j:j+8] for j in range(0,n,8)))):
                dest[8*i:8*i+len(packet)]=packet
            dest.flush();del dest;(args.output/'boundary.tmp').replace(disk)
        values=np.load(disk,mmap_mode='r',allow_pickle=False)
        used=set(cols);free=[j for j in range(width) if j not in used]
        assert len(free)==61
        # X is the fixed pivot block to the free columns. Augment the right hand
        # sides to compute the projection by one modular factorization.
        sq=nmod_mat(n,n,P); right=nmod_mat(n,len(free),P)
        for i,line in enumerate(values):
            for j,v in enumerate(line[cols]):sq[i,j]=int(v)
            for j,v in enumerate(line[free]):right[i,j]=int(v)
        det=int(sq.det());assert det==cert['determinant']
        solved=sq.solve(right)
        transition=np.asarray(solved.tolist(),dtype=np.int64)
        def residual(v):
            return (v[:,free]-matmul_mod(v[:,cols],transition,P))%P
        for start in range(0,n,144):assert not residual(values[start:start+144].astype(np.int64)).any()
        b=residual(rhs[None,:]); before=int(np.count_nonzero(b))
        extra=[];images=[];sp=RowSpace(len(free),P);tested=0
        for start in range(5000,5000+args.batches,32):
            if not b.any():break
            packets=list(pool.map(maps.sample,range(start,min(start+32,5000+args.batches))))
            ww=[w for ws,_ in packets for w in ws]; vv=np.vstack([v for _,v in packets]); tested+=len(ww)
            rem=residual(vv); keep=sp.extend(rem)
            extra.extend(ww[i] for i in keep);images.extend(vv[i].astype(np.uint32) for i in keep)
            assert n+len(extra)<=7254
            if sp.pivots and not ((b-b[:,sp.pivots]@sp.rows)%P).any(): break
        if extra:
            extra_mat=np.stack(images).astype(np.int64); raw=residual(extra_mat)
            h=nmod_mat(raw[:,sp.pivots].T.tolist(),P).solve(nmod_mat(b[:,sp.pivots].T.tolist(),P))
            coeff=np.array(h.tolist(),dtype=np.int64).ravel()
        else: extra_mat=np.zeros((0,width),dtype=np.int64);coeff=np.empty(0,dtype=np.int64)
        rem_rhs=(rhs-coeff@extra_mat)%P
        unresolved=residual(rem_rhs[None,:]);ok=not unresolved.any()
        record=dict(prime=P,target_profile=16,width=width,profile_dimension=7254,initial_minor_order=n,
            initial_minor_determinant=det,initial_residual_nonzero=before,extra_independent=len(extra),extra_tested=tested,
            augmented_rank_lower_bound=n+len(extra),rhs_terms=87234,rhs_sha256=sha(args.rhs),rhs_values_sha256=sha(args.output/'rhs.npy'),
            after_residual_nonzero=int(np.count_nonzero(unresolved)),projected_equation_solved=ok)
        if ok:
            x=sq.transpose().solve(nmod_mat(n,1,list(map(int,rem_rhs[cols])),P))
            a=np.asarray(x.tolist(),dtype=np.int64).ravel()
            replay=coeff@extra_mat%P
            for j in range(0,n,128):
                replay=(replay+a[j:j+128]@values[j:j+128].astype(np.int64))%P
            assert np.array_equal(replay,rhs)
            certificate=dict(base_minor_sha256=MINOR_HASH,base_coefficients=a.tolist(),extra_rows=extra,
                             extra_coefficients=coeff.tolist(),prime=P,profile=16)
            file=args.output/'solution.json.gz';file.write_bytes(gzip.compress(json.dumps(certificate,separators=(',',':')).encode(),mtime=0))
            record.update(replayed_coordinates=width,residual_nonzero=0,solution_sha256=sha(file),solution_bytes=file.stat().st_size,
                          base_nonzero=int(np.count_nonzero(a)),extra_nonzero=int(np.count_nonzero(coeff)))
        maps.write_json(args.output/'profile_result.json',record);print(json.dumps(record),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runtime',type=Path,default=maps.RUNTIME);p.add_argument('--plan',type=Path,
        default=maps.HERE/'classical_order1.json');p.add_argument('--minor',type=Path,default=maps.HERE/'other_maps/task_3/minor.json.gz')
    p.add_argument('--rhs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',choices=('pack','rhs','solve'),required=True)
    p.add_argument('--batches',type=int,default=1024);p.add_argument('--workers',type=int,default=16);a=p.parse_args();sys.path.insert(0,str(a.runtime.resolve()))
    if a.mode=='pack':pack_rhs(a.rhs,a.output)
    else:
        assert sha(a.plan)==PLAN_HASH and sha(a.minor)==MINOR_HASH and sha(a.rhs)==PROFILE_HASH
        cert=json.loads(gzip.decompress(a.minor.read_bytes()));plan=json.loads(a.plan.read_text())['incoming']
        assert cert['map']['targets']==[16] and cert['map']['prime']==P
        setup(a.runtime.resolve(),plan,[16],cert['map']['sources'],P);a.output.mkdir(parents=True,exist_ok=True)
        if a.mode=='rhs':rhs_values(a)
        else:solve(a,cert)
