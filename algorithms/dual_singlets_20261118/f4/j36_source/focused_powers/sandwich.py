"""Classical boundaries from traces with two separated open colour insertions."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import combinations,product
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'global_boundaries'))
from global_source_boundaries import minor, colour
from marked_jet_blocks import MarkEvaluation, branch_splits, group_degrees, derivative
from mixed_source_search import previous, solve_left, square_determinant
from adjoint_edge_search import central, old
from global_source_boundaries import HEIGHTS,subsets
from three_axis import Point, FULL
import three_axis

@cache
def splits3(mask,degree):
    available=[[i for i in range(3) if mask & (1<<(3*j+i))] for j in range(5)]
    return [(mask^sum(1<<(3*g+i) for g,s in enumerate(part) for i in s),
             sum(1<<(3*g+i) for g,s in enumerate(part) for i in s),
             (-1)**sum(i>j for a,s in zip(available,part) for i in a if i not in s for j in s))
            for part in product(*(combinations(a,k) for a,k in zip(available,degree)))]

HEIGHT={2:tuple(HEIGHTS),3:(3,)*5}
CUTS={
 2:((4,0,3,0),(4,4,0,0),(3,0,3,0),(4,0,2,0),(4,1,3,0),(4,3,0,0)),
 3:((3,3,0,0,0),(3,3,0,0,0),(3,3,1,0,0),(3,3,2,0,0),(3,2,1,0,0),(3,2,2,0,0))}
MODE=(('ad','jo'),('ad','ad'),('jo','jo'))
N=14
WIDTH=110

def make(rng,counts,marks):
    words=[]
    for n in range(N):
        groups=list(map(int,rng.permutation([j for j,v in enumerate(counts) for _ in range(v)])))
        sizes=[1]*len(groups) if marks==0 else [2]+[1]*(len(groups)-2) if marks==1 else (
            [3]+[1]*(len(groups)-3) if n%2==0 else [2,2]+[1]*(len(groups)-4))
        tokens=[];k=0
        for j in rng.permutation(sizes):
            j=int(j);tokens.append((('x','y','z')[j-1],*groups[k:k+j]));k+=j
        cut=int(rng.integers(1,len(tokens)))
        words.append((tuple(tokens[:cut]),tuple(tokens[cut:])))
    return words

@cache
def config(stage,order,case):
    cuts=CUTS[stage][case%6];other=tuple(a-b for a,b in zip(HEIGHT[stage],cuts))
    rng=np.random.default_rng(3809911+stage*1024+order*137+case)
    two=(order==2 and case%3!=1)
    return cuts,make(rng,cuts,1 if two else order),make(rng,other,int(two)),two,MODE[(case//6)%len(MODE)]

class Sandwich:
    def __init__(self,point,p,stage):
        self.e=point;self.p=p;self.stage=stage

    def evaluate(self,order,case):
        point=self.e;p=self.p;stage=self.stage
        cut,a,b,two,(f,g)=config(stage,order,case)
        pairs=(subsets(cut) if stage==2 else [(m,FULL^m,s) for m,s in three_axis.subsets(cut)])
        read=point.words() if stage==2 else point.values()
        first=point.frames[f][0].reshape(-1,26,26)
        trace_second=point.frames[g][1]
        z=np.zeros((len(first),trace_second.shape[1]))
        def wd(w):
            return tuple(sum(t[1:].count(k) for t in w) for k in range(len(cut)))
        @cache
        def piece(word,mask):
            left,right=word
            splits=(branch_splits if stage==2 else splits3)(mask,wd(right))
            result=z.copy()
            for l,r,s in splits:
                ll=read(left,l);rr=read(right,r)
                corners=(ll@first%p)@rr%p
                result+=s*(corners.reshape(len(first),676)@trace_second%p)
            return result%p
        diff=derivative if stage==2 else three_axis.derivative
        def d(spec,mask):
            l,r=spec
            ps=(-1)**sum(t[0]=='y' for t in l)
            return (sum((s*piece((v,r),mask) for s,v in diff(l)),z)
                    +ps*sum((s*piece((l,v),mask) for s,v in diff(r)),z))%p
        def contract(x,y):
            joined=(point.frames[f][2]@x%p)@point.frames[g][2]%p
            return joined.reshape(len(x),-1)@y.reshape(len(y),-1).T%p
        total=np.zeros((len(a),len(b)))
        for m,n,sg in pairs:
            u=np.stack([d(w,m) for w in a]);v=np.stack([piece(w,n) for w in b])
            val=contract(u,v)
            if two:
                val-=contract(np.stack([piece(w,m) for w in a]),np.stack([d(w,n) for w in b]))
            total+=sg*val
            piece.cache_clear();read.cache_clear()
        return np.asarray(total.ravel()%p,dtype=np.uint16)

@cache
def indices(stage,order,c):
    p=1009
    with np.errstate(all='ignore'):
        point=MarkEvaluation(16654,p) if stage==2 else Point(88422,p)
        row=Sandwich(point,p,stage).evaluate(order,c)
    live=np.flatnonzero(row)
    return np.sort(np.random.default_rng(20261124+143*stage+15*order+c).choice(live,min(len(live),WIDTH),replace=False))

@cache
def first(stage,order):
    if stage==2:return previous()['rows']+central()['new_rows'],3002
    z=json.loads((HERE/('three-1-1009-429.json' if order==1 else 'three-2-1009-1664.json')).read_bytes())
    return z['rows'],z['rank']

def run_old(args):
    stage,order,i,p=args
    if stage==2:return old((i,p))
    nm='three-1-1009-429.json' if order==1 else 'three-2-1009-1664.json'
    z=json.loads((HERE/nm).read_bytes())
    return three_axis.point((i,p,order,range(len(three_axis.CUTS))))[z['columns']]

def run_new(args):
    stage,order,i,p,case=args
    with np.errstate(all='ignore'):
        point=MarkEvaluation(14337+i,p) if stage==2 else Point(i,p)
        return Sandwich(point,p,stage).evaluate(order,case)[indices(stage,order,case)]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['search','controls','replay'])
    ap.add_argument('--stage',type=int,choices=[2,3],required=True);ap.add_argument('--order',type=int,choices=[1,2],required=True)
    ap.add_argument('--prime',type=int,default=1009,choices=[1009,1013]);ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--cases',type=int,default=18);ap.add_argument('--certificate',type=Path);args=ap.parse_args()
    stage=args.stage;order=args.order;p=args.prime;assert stage!=2 or order==2
    if args.mode=='controls':
        with np.errstate(all='ignore'):
            pt=MarkEvaluation(95611,p,True) if stage==2 else Point(81,p,True)
            assert all(not Sandwich(pt,p,stage).evaluate(2,c).any() for c in range(0,args.cases))
        print('sandwich closure controls passed',stage,p,flush=True);return
    orig,n=first(stage,order)
    dest=Path(f'sandwich-{stage}-{order}.json')
    saved=json.loads(args.certificate.read_bytes()) if args.mode=='replay' else None
    ext=saved['new_rows'] if saved else list(range(7300,7300+(32 if stage==2 else 65)))
    ids=orig+ext
    assert len(set(ids))==len(ids)
    maximum=3017 if stage==2 else (428 if order==1 else 1862-419)
    with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
        base=np.array(list(pool.map(run_old,((stage,order,i,p) for i in ids),chunksize=2)))
        with np.errstate(all='ignore'):w,det=solve_left(base[:n],base[n:],p)
        baseline={(2,2,1009):833,(2,2,1013):615,(3,1,1009):200,(3,1,1013):723,(3,2,1009):680,(3,2,1013):780}
        assert det==baseline[stage,order,p]
        print('base',stage,order,n,det,flush=True)
        keep=[]; residual=np.zeros((len(ext),0))
        if saved:
            keep=saved['new_columns'];needed=sorted({i//WIDTH for i in keep})
            res=np.empty((len(ext),len(keep)))
        else:
            needed=[c for k in range(6) for c in (k,6+k,12+k) if c<args.cases]
        for case in needed:
            if not len(indices(stage,order,case)):continue
            raw=np.array(list(pool.map(run_new,((stage,order,i,p,case) for i in ids),chunksize=2)))
            with np.errstate(all='ignore'):
                reduction=(raw[n:]-w@raw[:n])%p
                if saved:
                    for col,j in enumerate(keep):
                        if j//WIDTH==case:res[:,col]=reduction[:,j%WIDTH]
                    continue
                union=np.column_stack([residual,reduction]);rr,cc,dd=minor(union,p)
            cols=keep+[case*WIDTH+j for j in range(reduction.shape[1])]
            keep=[cols[j] for j in cc];residual=union[:,cc]
            assert n+len(rr)<=maximum
            record=dict(stage=stage,order=order,prime=p,previous_rank=n,rank=n+len(rr),base_determinant=det,
                        complement_determinant=int(dd),determinant=int(det*dd%p),new_rows=[ext[j] for j in rr],
                        new_columns=keep,last_family=case)
            temp=dest.with_suffix('.part');temp.write_text(json.dumps(record)+'\n');temp.replace(dest)
            print('family',case,record['rank'],record['determinant'],flush=True)
            if record['rank']==maximum:break
        if saved:
            with np.errstate(all='ignore'):ds=square_determinant(res,p)
            record=dict(stage=stage,order=order,prime=p,rank=n+len(keep),base_determinant=det,complement_determinant=int(ds),determinant=int(ds*det%p))
            dest.with_name(f'sandwich-replay-{stage}-{order}.json').write_text(json.dumps(record)+'\n')
            print(record,flush=True)

if __name__=='__main__':main()
