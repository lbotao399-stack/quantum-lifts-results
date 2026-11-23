"""Test the remaining classical source with two adjoint colour edges."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from pathlib import Path
import json
import numpy as np
from multiprocessing import get_context

from global_source_boundaries import colour,HEIGHTS,HERE,minor,write,subsets
from marked_jet_blocks import marked_extensions
from mixed_colour_edges import EdgeEvaluation,differentiated
from mixed_source_search import (bank,base_point,new_point,previous,solve_left,square_determinant)

N=38
WIDTH=480
CUTS=((4,0,3,0),(4,4,0,0),(3,0,3,0),(4,0,2,0),(4,3,0,0),(4,1,3,0),
      (4,0,2,1),(4,1,2,0),(3,0,3,1),(3,1,3,0),(4,2,0,1),(4,2,1,0))


@cache
def bridge(p):
    generators,_=colour(p)
    return (generators[None]@generators[:,None]%p).transpose(0,1,3,2).reshape(2704,676).T.astype(float)


class AdjointEdges(EdgeEvaluation):
    def __init__(self,seed,p,brackets=False):
        super().__init__(seed,p,brackets)
        self.adbridge=bridge(p);self.squarezero=np.zeros((52,52))

    def build(self,variant):
        direct=self.words()
        @cache
        def adj_prefix(w,mask):
            if not w:return np.eye(52)
            terms=marked_extensions(w[-1],mask)
            if not terms:return self.squarezero
            left=np.concatenate([adj_prefix(w[:-1],m)*sg for m,sg,k in terms],axis=1)
            right=np.concatenate([self.maps[k] for m,sg,k in terms])
            return left@right%self.p
        @cache
        def process(w,mask):
            if w[-1]=='cross':
                if variant==0:
                    return (direct(w[:-1],mask).reshape(676)@self.adbridge%self.p).reshape(52,52)
                # The adjoint representation supplies an additional square colour word.
                a=adj_prefix(w[:-1],mask)
                # The invariant trace metric identifies the first open adjoint leg with its dual.
                return self.trace_metric@a%self.p
            side,role,*groups=w[-1]
            terms=marked_extensions((role,*groups),mask)
            if not terms:return self.squarezero
            if side=='adj':
                left=np.concatenate([process(w[:-1],m)*sg for m,sg,k in terms],axis=1)
                right=np.concatenate([self.maps[k] for m,sg,k in terms])
            else:
                left=np.concatenate([self.maps[k].T*sg for m,sg,k in terms],axis=1)
                right=np.concatenate([process(w[:-1],m) for m,sg,k in terms])
            return left@right%self.p
        def clear():
            direct.cache_clear();process.cache_clear();adj_prefix.cache_clear()
        return process,clear

    def evaluate(self,c,indices):
        cuts,a,b,external,variants=spec(c)
        i,j=np.unravel_index(indices,(N,N))
        sa,ia=np.unique(i,return_inverse=True);sb,ib=np.unique(j,return_inverse=True)
        a=[a[z] for z in sa];b=[b[z] for z in sb]
        ra,ca=self.build(variants[0]);rb,cb=self.build(variants[1])
        def diff(words,mask,reader):
            return np.stack([sum((s*reader(w,mask) for s,w in differentiated(v)),self.squarezero)%self.p for v in words])
        gram=self.frames['ad'][2]
        def pairs(x,y):
            u=(gram@x%self.p)@gram%self.p
            return (u.reshape(len(x),2704)@y.reshape(len(y),2704).T)[ia,ib]%self.p
        total=np.zeros(len(indices))
        for m,n,sign in subsets(cuts):
            total+=sign*pairs(diff(a,m,ra),np.stack([rb(w,n) for w in b]))
            if external:total-=sign*pairs(np.stack([ra(w,m) for w in a]),diff(b,n,rb))
            ca();cb()
        return np.asarray(total%self.p,dtype=np.uint16)

    @property
    def trace_metric(self):
        return trace_metric(self.p)


@cache
def trace_metric(p):
    g,_=colour(p);return (g.reshape(52,676)@g.transpose(0,2,1).reshape(52,676).T%p).astype(float)


@cache
def spec(c):
    n=len(CUTS);cuts=CUTS[c%n];other=tuple(h-k for h,k in zip(HEIGHTS,cuts))
    rng=np.random.default_rng(7811367+c);external=c%3!=1
    return cuts,bank(rng,cuts,N,1 if external else 2),bank(rng,other,N,1 if external else 0),external,((c//n)%2,0)


@cache
def screens():
    return [np.asarray(v,dtype=int) for v in json.loads((HERE/'adjoint_screen.json').read_bytes())['indices']]


@cache
def central():
    return json.loads((HERE/'mixed-minor.json').read_bytes())


def old(task):
    i,p=task
    return np.r_[base_point(task),new_point((i,p,central()['new_columns']))]


def new(task):
    i,p,cols=task
    with np.errstate(all='ignore'):
        e=AdjointEdges(14337+i,p)
        if cols is None:
            return np.concatenate([e.evaluate(c,v) for c,v in enumerate(screens())])
        cols=np.asarray(cols,dtype=int);out=np.empty(len(cols),dtype=np.uint16)
        for c,indices in enumerate(screens()):
            take=np.flatnonzero(cols//WIDTH==c)
            if len(take):out[take]=e.evaluate(c,indices[cols[take]%WIDTH])
        return out


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=('screen','controls','search','replay'))
    ap.add_argument('--workers',type=int,default=1);ap.add_argument('--prime',type=int,choices=(1009,1013),default=1009)
    ap.add_argument('--certificate',type=Path,default=HERE/'adjoint-minor.json');args=ap.parse_args()
    assert __debug__ and args.workers>=1;p=args.prime
    if args.action=='screen':
        data=[]
        with np.errstate(all='ignore'):
            e=AdjointEdges(714222,1009)
            for c in range(2*len(CUTS)):
                vals=e.evaluate(c,np.arange(N*N));live=np.flatnonzero(vals);rng=np.random.default_rng(771029+c)
                pick=rng.choice(live,WIDTH,replace=False) if len(live)>=WIDTH else np.r_[live,rng.choice(np.flatnonzero(vals==0),WIDTH-len(live),replace=False)]
                data.append(sorted(map(int,pick)))
        write(HERE/'adjoint_screen.json',dict(seed=714222,indices=data))
        return
    if args.action=='controls':
        with np.errstate(all='ignore'):
            e=AdjointEdges(15690,p,brackets=True)
            for c,v in enumerate(screens()): assert not e.evaluate(c,v).any()
        print('All two-adjoint boundary controls passed.',flush=True);return
    save=json.loads(args.certificate.read_bytes()) if args.action=='replay' else None
    extras=save['new_rows'] if save else list(range(7000,7040))
    ids=previous()['rows']+central()['new_rows']+extras
    assert len(set(ids))==len(ids) and central()['rank']==3002
    cols=save['new_columns'] if save else None
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        base=np.asarray(list(pool.map(old,((i,p) for i in ids),chunksize=2)))
        with np.errstate(all='ignore'):weights,d=solve_left(base[:3002],base[3002:],p)
        if p==1009:assert d==833
        print('retained',3002,d,flush=True)
        # A zero-size residual determinant certifies the previous rank only.
        if cols==[]:
            print(json.dumps(dict(prime=p,size=3002,base_minor=d,complement_minor=1,determinant=d)),flush=True);return
        raw=np.asarray(list(pool.map(new,((i,p,cols) for i in ids),chunksize=2)),dtype=float)
    with np.errstate(all='ignore'):small=(raw[3002:]-weights@raw[:3002])%p
    if save:
        ds=square_determinant(small,p)
        print(json.dumps(dict(prime=p,size=len(ids),base_minor=d,complement_minor=ds,determinant=int(d*ds%p))),flush=True)
    else:
        r,c,ds=minor(small,p)
        assert 3002+len(r)<=3017
        write(args.certificate,dict(prime=p,previous_rank=3002,new_rows=[extras[i] for i in r],new_columns=c,
                                    complement_minor=ds,rank=3002+len(r),determinant=int(d*ds%p)))
        print('extended',3002+len(r),d*ds%p,flush=True)


if __name__=='__main__':main()
