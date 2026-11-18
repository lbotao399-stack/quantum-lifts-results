"""SO(8) BMN: trace graviton products and polarized Pfaffian descendants.

Coefficients and Qplus are unchanged from the pinned upstream convention.
Pfaffian expressions use the existing raw epsilon contraction (384 Pf).
Products are evaluated directly, so two epsilon factors need no expansion.
This tests a specified closed span; a positive residual is only a candidate.
"""
import os,sys,pathlib,json,time,argparse,signal,resource
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
ROOT=pathlib.Path(__file__).resolve().parent
CORE=ROOT.parent/'so8-bmn'
sys.path[:0]=[str(CORE/'src'),str(CORE/'vendor'),str(CORE/'deps')]
import numpy as np
from scipy.sparse import csr_matrix
from collections import defaultdict
from itertools import combinations,combinations_with_replacement
from math import comb,ceil
from pilot import expected,manifest
from invariants import canon,count,sector_basis,action,action_expr,clear_caches
from evaluate import Frame
from finite_field import rank_mod
from fortuity.bmn.words import graviton_generators,charges_of_content


def graph_expr(expr):
    out=defaultdict(int)
    for traces,c in expr.items():
        s,w=canon((),traces)
        if s:out[w]+=s*c
    return {w:c for w,c in out.items() if c}


def expr_sec(expr):
    secs={charges_of_content(count(w)) for w in expr}
    assert len(secs)==1,secs
    return secs.pop()


def generators(jmax):
    out=[];seen=set()
    def add(name,e,reflection):
        if not e:return
        sec=expr_sec(e)
        if 2*sum(sec[1:])>jmax:return
        key=tuple(sorted(e.items()))
        if key in seen:return
        seen.add(key)
        degrees={sum(count(w)[3:6])%2 for w in e}
        assert len(degrees)==1
        out.append({'name':name,'expr':e,'sector':sec,'odd':degrees.pop(),'reflection':reflection})
    for name,expr in graviton_generators(jmax//2,1):
        add(name,graph_expr(expr),0)
    for inds in combinations_with_replacement(range(3),4):
        s,w=canon(tuple((i,) for i in inds),())
        e={w:s}
        for n in range(4):
            for actions in combinations(range(3),n):
                d=e
                for a in actions:d=action_expr(d,str(a))
                add(('Pf',inds,actions),d,1)
    return out


def monomials(gens,sec,reflection):
    ids=[i for i,g in enumerate(gens) if all(x<=y for x,y in zip(g['sector'],sec))]
    out=[]
    def rec(start,rem,par,chosen):
        if not any(rem):
            if par==reflection:out.append(tuple(chosen))
            return
        for pos in range(start,len(ids)):
            i=ids[pos];g=gens[i]
            if all(x<=y for x,y in zip(g['sector'],rem)):
                rec(pos+g['odd'],tuple(y-x for x,y in zip(g['sector'],rem)),
                    par^g['reflection'],chosen+[i])
    rec(0,tuple(sec),0,[])
    return out


def eval_expr(expr,fr,maxdegree):
    out={}
    for w,c in expr.items():
        d=sum(count(w)[3:6])
        if d>maxdegree:continue
        value=fr.invariant(w).reshape(fr.P,comb(fr.k,d))*(c%fr.p)%fr.p
        out[d]=(out[d]+value)%fr.p if d in out else value
    return out


def eval_products(gens,monos,fr,degree):
    cache={i:eval_expr(gens[i]['expr'],fr,degree) for m in monos for i in m}
    rows=np.zeros((len(monos),fr.ncols),dtype=np.int64)
    for row,m in enumerate(monos):
        value={0:np.ones((fr.P,1),dtype=np.int64)}
        for i in m:
            nxt={}
            for d,v in value.items():
                for e,w in cache[i].items():
                    if d+e>degree:continue
                    product,k=fr.scal_mul(v,d,w,e)
                    nxt[k]=(nxt[k]+product)%fr.p if k in nxt else product
            value=nxt
        if degree in value:rows[row]=value[degree].reshape(-1)
    return rows


def check_generators(gens):
    # {Q,Qplus_m} is the infinitesimal gauge action of phi_m.
    # For epsilon graphs gauge identities are not part of graph canonization;
    # verify the resulting polynomial identities in generic skew matrices.
    tested=0
    for p in (1000003,1000033):
        for k,g in enumerate(gens):
            image=action_expr(g['expr'])
            for d in {sum(count(w)[3:6]) for w in image}:
                fr=Frame(d,2,9817+31*k,p,k=d+2)
                value=eval_expr(image,fr,d).get(d)
                assert value is None or not value.any(),('generator not closed',g['name'],p,d)
            tested+=1
    return {'generator_evaluation_checks':tested,'primes':[1000003,1000033]}


def sector(sec,parity,gens,known,prime=1000003,seed=1701):
    started=time.monotonic();exact=expected(sec,parity);dim=sum(exact.values())
    rank_in=known[(tuple(sec),(-1)**parity)]['rank_in']
    outgoing=(sec[0]-1,)+tuple(sec[1:])
    rank_out=known.get((outgoing,(-1)**parity),{}).get('rank_in',0)
    h=dim-rank_in-rank_out
    result={'sector':sec,'reflection_parity':(-1)**parity,'prime':prime,'seed':seed,
            'dimension':dim,'rank_in':rank_in,'rank_out':rank_out,'cohomology':h}
    assert h>=0
    if not h:
        return dict(result,protected=0,residual=0,seconds=time.monotonic()-started,
                    method='cohomology vanishes by classical rank saturation')
    monos=monomials(gens,sec,parity)
    words=sector_basis(sec,parity);incoming=sector_basis((sec[0]+1,)+tuple(sec[1:]),parity)
    idx={w:i for i,w in enumerate(words)};degree=np.array([sum(count(w)[3:6]) for w in words])
    rr=[];cc=[];vv=[]
    for i,w in enumerate(incoming):
        for v,c in action(w).items():rr.append(i);cc.append(idx[v]);vv.append(c)
    Q=csr_matrix((np.asarray(vv,dtype=np.int64),(rr,cc)),shape=(len(incoming),len(words)))
    boundaries=[];protected=[];evaluation={}
    for d,m in sorted(exact.items()):
        if not m:continue
        ids=np.flatnonzero(degree==d);selected=[words[i] for i in ids]
        blocks=[];gblocks=[];used=0;k=0 if d==0 else d+2;perpoint=comb(k,d)
        for attempt in range(4):
            target=m+8 if attempt==0 else 2*used
            while used<target:
                npts=min(max(1,32//perpoint),8,ceil((target-used)/perpoint))
                fr=Frame(d,npts,seed+100003*d+7919*used,prime,k=k)
                blocks.append(fr.values(selected).astype(np.int32))
                gblocks.append(eval_products(gens,monos,fr,d).astype(np.int32));used+=fr.ncols
                del fr
            X=np.concatenate(blocks,axis=1).astype(np.int64)
            r=rank_mod(X,prime)
            if r==m:break
            assert r<m
        assert r==m,('incomplete target evaluation',sec,parity,d,r,m)
        boundaries.append(np.asarray(Q[:,ids]@X)%prime)
        protected.append(np.concatenate(gblocks,axis=1).astype(np.int64))
        evaluation[d]={'exact_dimension':m,'evaluation_rank':r,'columns':used}
    B=np.concatenate(boundaries,axis=1);G=np.concatenate(protected,axis=1)
    rin=rank_mod(B,prime);assert rin==rank_in,(sec,parity,rin,rank_in)
    combined=rank_mod(np.concatenate([B,G],axis=0),prime)
    protected_rank=combined-rin;residual=h-protected_rank
    assert residual>=0,('closed span exceeds cohomology',sec,parity,h,protected_rank)
    return dict(result,protected=protected_rank,residual=residual,combined_rank=combined,
                generator_products=len(monos),evaluation=evaluation,seconds=time.monotonic()-started,
                maxrss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                method='faithful full cochain evaluation; rank(boundaries + trace/Pfaffian products)')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--j',type=int,required=True)
    ap.add_argument('--prime',type=int,default=1000003);ap.add_argument('--seed',type=int,default=1701)
    ap.add_argument('--seconds',type=int,default=1200);ap.add_argument('--output',required=True)
    ap.add_argument('--check',action='store_true');a=ap.parse_args()
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('protected span wall limit')))
    signal.alarm(a.seconds)
    knownrows=[json.loads(s) for s in (CORE/'runs_flint'/('j%d_p%d.jsonl'%(a.j,a.prime))).read_text().splitlines()]
    known={(tuple(r['sector']),r['reflection_parity']):r for r in knownrows}
    gens=generators(a.j)
    if a.check:print(json.dumps(check_generators(gens)),flush=True)
    p=pathlib.Path(a.output);p.parent.mkdir(parents=True,exist_ok=True)
    complete=set()
    if p.exists():complete={(tuple(r['sector']),r['reflection_parity']) for r in map(json.loads,p.read_text().splitlines())}
    with p.open('a') as f:
        for sec,par in manifest(a.j):
            if (tuple(sec),(-1)**par) in complete:continue
            r=sector(sec,par,gens,known,a.prime,a.seed)
            f.write(json.dumps(r)+'\n');f.flush()
            print(json.dumps({k:r[k] for k in ('sector','reflection_parity','cohomology','protected','residual','seconds')}),flush=True)
            clear_caches()


if __name__=='__main__':main()
