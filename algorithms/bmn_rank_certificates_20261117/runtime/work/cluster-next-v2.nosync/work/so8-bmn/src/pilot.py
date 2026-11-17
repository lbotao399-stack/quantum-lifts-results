#!/usr/bin/env python3
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import sys,pathlib,time,json,argparse,resource,signal
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'vendor'))
import numpy as np
from collections import defaultdict
from scipy.sparse import csr_matrix
from fortuity.bmn.words import contents_of_sector
from finite_field import rank_mod
from math import comb,ceil
from exact_dims import multiplicity,scalar_ring_check
from invariants import sector_basis,count,action,action_expr
from evaluate import Frame


def expected(sec,parity):
    bydegree=defaultdict(int)
    for c in contents_of_sector(sec):
        bydegree[sum(c[3:6])]+=multiplicity(c,parity)
    return dict(bydegree)


def sector(sec,parity,prime=1000003,seed=1701):
    started=time.monotonic();times={}
    exact=expected(sec,parity)
    m=sum(exact.values());times['dimensions']=time.monotonic()-started
    base={'sector':sec,'reflection_parity':1 if parity==0 else -1,'prime':prime,'seed':seed,
          'rank_backend':'python-flint nmod_mat',
          'dimension':m,'dimensions_by_fermions':exact}
    if not m:
        return dict(base,rank_in=0,basis_words=0,spanning_verified=True,seconds=time.monotonic()-started,timing=times)
    t=time.monotonic()
    words=sector_basis(sec,parity)
    incoming=sector_basis((sec[0]+1,)+tuple(sec[1:]),parity)
    idx={w:i for i,w in enumerate(words)}
    degree=np.array([sum(count(w)[3:6]) for w in words])
    times['enumerate']=time.monotonic()-t;t=time.monotonic()
    rr=[];cc=[];vv=[]
    for i,w in enumerate(incoming):
        for v,c in action(w).items():
            assert v in idx,('missing target',sec,v)
            rr.append(i);cc.append(idx[v]);vv.append(c)
    Q=csr_matrix((np.array(vv,dtype=np.int64),(rr,cc)),shape=(len(incoming),len(words)))
    times['Q']=time.monotonic()-t;t=time.monotonic()
    # Each target degree has a separate Grassmann frame. Q-images can mix
    # degrees: concatenate their evaluations before taking one total rank.
    images=[];ranks={};points={};cols={};generators={}
    for d,dim in sorted(exact.items()):
        if not dim:continue
        ids=np.flatnonzero(degree==d);selected=[words[i] for i in ids]
        blocks=[];used=0;npoints=0
        k=0 if d==0 else d+2
        per_point=comb(k,d)
        for attempt in range(4):
            target=dim+8 if attempt==0 else 2*used
            while used<target:
                npts=min(max(1,32//per_point),8,ceil((target-used)/per_point))
                fr=Frame(d,npts,seed+100003*d+7919*used,prime,k=k)
                blocks.append(fr.values(selected).astype(np.int32));used+=fr.ncols;npoints+=npts
                del fr
            X=np.concatenate(blocks,axis=1).astype(np.int64)
            rank=rank_mod(X,prime)
            if rank==dim:break
            assert rank<dim,('dimension mismatch',sec,parity,d,rank,dim)
        assert rank==dim,('incomplete evaluation',sec,parity,d,rank,dim)
        ranks[d]=rank;points[d]=npoints;cols[d]=used;generators[d]=k
        images.append(np.asarray(Q[:,ids]@X)%prime)
    B=np.concatenate(images,axis=1) if images else np.zeros((len(incoming),0),dtype=np.int64)
    times['evaluate']=time.monotonic()-t;t=time.monotonic()
    rank_in=rank_mod(B,prime)
    times['rank']=time.monotonic()-t
    return dict(base,rank_in=rank_in,basis_words=len(words),incoming_words=len(incoming),
                spanning_verified=True,rank_by_fermions=ranks,evaluation_points=points,
                evaluation_columns=cols,grassmann_generators=generators,
                seconds=time.monotonic()-started,timing=times,
                maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                certification='exact integer cochain dimension; faithful evaluation and ranks over specified finite field')


def columns(j):
    total=j//2
    return [(i,k,total-i-k) for i in range(total+1) for k in range(i+1)
            if 0<=total-i-k<=k]


def manifest(j):
    out=[]
    for q in columns(j):
        for n in range(j//3+1):
            sec=(n,)+q
            if contents_of_sector(sec):
                for parity in (0,1):out.append((sec,parity))
    return out


def checks():
    from invariants import canon,trace_canon,path_canon,invariants_content
    scalar=scalar_ring_check(12)
    checked=0
    for j in range(0,13,2):
        for sec,parity in manifest(j):
            for w in sector_basis(sec,parity)[:20]:
                assert not action_expr(action(w)),('Q squared',w)
                checked+=1
    # Direct normalization and reflection controls of the Pfaffian evaluator.
    fr=Frame(0,2,321,1000003)
    diagonal=np.zeros((2,1,8,8),dtype=np.int64)
    for i,v in enumerate((2,3,5,7)):
        diagonal[:,:,2*i,2*i+1]=v;diagonal[:,:,2*i+1,2*i]=-v
    fr.L[0]=(0,diagonal%fr.p)
    pf=fr.invariant((((0,),)*4,()))
    assert np.all(pf==(384*2*3*5*7)%fr.p),pf
    assert not trace_canon((0,0,0))[0]
    assert not path_canon((0,0))[0]
    # Small gauge groups are independent scalar-ring checks of the character method.
    for n in range(9):
        assert multiplicity((n,0,0,0,0,0,0),0,3)==int(n%2==0)
        assert multiplicity((n,0,0,0,0,0,0),1,3)==0
    return {'scalar_ring_degree_12':scalar,'Q_squared_invariants_checked':checked,
            'pfaffian_normalization':384*2*3*5*7,'SO3_scalar_ring':True}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true')
    ap.add_argument('--sector');ap.add_argument('--parity',type=int,default=0)
    ap.add_argument('--prime',type=int,default=1000003);ap.add_argument('--output')
    ap.add_argument('--seconds',type=int,default=115)
    a=ap.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('wall time limit')))
    signal.alarm(a.seconds)
    if a.check:result=checks()
    else:result=sector(tuple(map(int,a.sector.split(','))),a.parity,a.prime)
    text=json.dumps(result,sort_keys=True)
    if a.output:pathlib.Path(a.output).write_text(text+'\n')
    print(text,flush=True)


if __name__=='__main__':main()
