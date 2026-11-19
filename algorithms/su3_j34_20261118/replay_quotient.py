"""Replay the profile-10 quotient minor, original training equations, or fresh trace evaluation.

The first two checks use only NumPy and hashed NPY inputs. The fresh check
uses the saved trace recipes and field evaluator without any training arrays.
All three compare actual coordinates, not just the number of mismatches.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

P = 1000003
SEED = 635729813


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest() if hasattr(hashlib,'file_digest') else hashlib.sha256(stream.read()).hexdigest()


def hv(a):
    return hashlib.sha256(np.asarray(a,dtype='<u4').tobytes()).hexdigest()


def state(folder):
    assert __debug__
    rec=json.loads((folder/'result.json').read_bytes()); c=np.load(folder/'coefficients.npy', allow_pickle=False)
    assert sha(folder/'result.json')==RESULT_HASH and sha(folder/'coefficients.npy')==rec['coefficients_file_sha256']==COEFF_HASH
    assert c.shape==(7616,) and c.dtype.str=='<u4' and c.max()<P and rec['prime']==P
    return rec,c.astype(np.int64)


def factor(a):
    """Pivot a rectangular matrix column by column; return rank and signed pivot product."""
    u=np.array(a,dtype=np.int64,copy=True); rank=0; det=1
    for j in range(u.shape[1]):
        hit=np.flatnonzero(u[rank:,j]);
        if not hit.size: continue
        i=rank+int(hit[0]);
        if i!=rank:
            u[[rank,i]]=u[[i,rank]]; det=-det
        d=int(u[rank,j]); det=det*d%P
        m=u[rank+1:,j]*pow(d,-1,P)%P
        u[rank+1:,j:]=(u[rank+1:,j:]-m[:,None]*u[rank,j:])%P
        rank+=1
        if rank==len(u): break
    return rank,det%P


def local(folder):
    r,c=state(folder); q,b=(np.load(folder/n,allow_pickle=False) for n in ('quotient.npy','quotient_rhs.npy'))
    assert q.shape==(768,365) and b.shape==(365,) and q.dtype.str==b.dtype.str=='<u4' and max(q.max(),b.max())<P
    assert (hv(q),hv(b))==(r['quotient_data_sha256'],r['quotient_rhs_data_sha256'])
    rows,cols=r['quotient_pivot_rows'],r['quotient_pivot_cols'];assert len(set(rows))==len(rows)==365 and sorted(cols)==list(range(365))
    rank,det=factor(q[np.ix_(rows,cols)])
    assert (rank,det)==(r['quotient_rank'],r['quotient_determinant'])==(365,143392)
    assert sha(folder/'words.json.gz')==r['word_file_sha256'];labels=[s for s,_ in json.loads(gzip.decompress((folder/'words.json.gz').read_bytes()))];assert len(labels)==len(q)
    by_profile={s:factor(q[np.array(labels)==int(s)])[0] for s in r['source_stats']}
    assert all(by_profile[s]==r['source_stats'][s]['rank'] for s in by_profile)
    lhs=c[6848:]@q.astype(np.int64)%P; np.testing.assert_array_equal(lhs,b)
    wide=json.loads((folder.parent/'wide/wide.json').read_bytes()); assert sha(folder.parent/'wide/wide.json')==r['old_record_sha256']; combined=wide['minor_determinant']*det%P;assert combined==r['new_minor_determinant']==102708
    return dict(prime=P,quotient_rank=rank,quotient_minor_determinant=det,ranks_by_profile=by_profile,combined_minor_determinant_in_recorded_pivot_order=combined,coordinates_checked=len(b),nonzero_residual=0,coefficients_file_sha256=COEFF_HASH,matched_rhs_sha256=hv(lhs))


def training(folder,wide):
    r,c=state(folder); old=json.loads((wide/'wide.json').read_bytes()); assert sha(wide/'wide.json')==r['old_record_sha256']
    rhs=np.load(wide/'rhs.npy',allow_pickle=False); assert rhs.shape==(7000,) and rhs.dtype.str=='<u4' and rhs.max()<P and hv(rhs)==old['rhs_sha256']
    out=np.zeros(len(rhs),dtype=np.int64); off=0; inputs={}
    for path,num,want in ((wide/'matrix.npy',6848,old['matrix_sha256']),(folder/'new_rows.npy',768,r['new_rows_data_sha256'])):
        a=np.load(path,mmap_mode='r',allow_pickle=False); assert a.shape==(num,len(rhs)) and a.dtype.str=='<u4'
        h=hashlib.sha256();
        for i in range(0,num,29):
            block=a[i:i+29]; h.update(block.tobytes()); assert block.max()<P and len(block)*(P-1)**2+P < np.iinfo(np.int64).max
            out=(out+c[off+i:off+i+len(block)]@block.astype(np.int64))%P
        assert h.hexdigest()==want; off+=num; inputs[path.name]=want
    assert off==len(c);np.testing.assert_array_equal(out,rhs)
    return dict(prime=P,producer_slurm_job='5650093',rows_multiplied=off,nonzero_coefficients=int(np.count_nonzero(c)),coordinates_checked=len(rhs),nonzero_residual=0,matched_rhs_sha256=hv(out),coefficient_file_sha256=COEFF_HASH,input_array_data_sha256=inputs)


def fresh(root,folder):
    import profile10_validation as v
    sys.path.insert(0,str(v.maps.RUNTIME))
    import su_trace_adapter_v1 as alg
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate as derivative
    from native_su_values_v1 import evaluate
    r,c=state(folder); ref=json.loads((folder/'fresh.json').read_bytes());assert sha(folder/'fresh.json')==FRESH_HASH and ref['frozen_coefficients_sha256']==COEFF_HASH and ref['seed']==r['precommitted_holdout_seed']==SEED
    p,words,_=v.inputs(root,root/'allcoord');fr,charge=v.frame(p,12,SEED);target=v.maps.pattern(p['targets'][10]); f=folder.parent
    for path,want in ((f/'new_words.json.gz',json.loads((f/'build.json').read_bytes())['new_words_sha256']),(folder/'words.json.gz',r['word_file_sha256'])):
        assert sha(path)==want; words+= [(s,v.maps.decode(w)) for s,w in json.loads(gzip.decompress(path.read_bytes()))]
    assert len(words)==len(set(words))==len(c);lhs=np.zeros(fr.ncols,dtype=np.int64);rhs=lhs.copy(); terms=list(v.rhs_terms(root,p).items())
    for batch in v.batches(terms,137): rhs=(rhs+evaluate([dict(batch)],fr)[0])%P
    for s in v.EDGE:
        active=[i for i,(t,w) in enumerate(words) if s==t and c[i]]
        for ids in v.batches(active,13):
            ex=[project(alg.parity_project({words[i][1]:1},1),charge,P) for i in ids]; lhs=(lhs+c[ids]@derivative(ex,fr,v.maps.pattern(p['sources'][s]),target))%P
    d=(rhs-lhs)%P
    for n,x in [('lhs',lhs),('rhs',rhs),('residual',d)]:np.testing.assert_array_equal(x,ref[n])
    assert hv(d)==ref['residual_sha256'] and np.count_nonzero(d)==ref['nonzero']==len(d)==120
    return dict(prime=P,seed=SEED,coordinates_checked=len(d),nonzero_residual=int(np.count_nonzero(d)),trace_terms=len(terms),active_coefficients=int(np.count_nonzero(c)),coefficients_file_sha256=COEFF_HASH,lhs_sha256=hv(lhs),rhs_sha256=hv(rhs),residual_sha256=hv(d),first_rhs=int(rhs[0]),first_lhs=int(lhs[0]),first_residual=int(d[0]))


# Pinned to the terminal production receipt, before replay.
RESULT_HASH = 'ef37638c5392cf9f2aa70ed345d8365393a1fc17ddd6ba6e3b5033d95bbb354c'
COEFF_HASH = '381a1bf86bc55850b3ed6df64a3bfc47a0456de1cfa576ffd47ac3d2df8b4dfb'
FRESH_HASH = '5a532ecb833b126a2a61b844291371cf1f18b93ef2d14353e768cecd482eeef1'

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['minor','training','fresh']);ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent);ap.add_argument('--quotient',type=Path);ap.add_argument('--wide',type=Path);ap.add_argument('--output',type=Path);args=ap.parse_args();q=args.quotient or args.root/'four_source/quotient'
    data=local(q) if args.mode=='minor' else training(q,args.wide or q.parent/'wide') if args.mode=='training' else fresh(args.root,q)
    text=json.dumps(data,indent=2)+'\n';print(text,end='');
    if args.output:args.output.write_text(text)
