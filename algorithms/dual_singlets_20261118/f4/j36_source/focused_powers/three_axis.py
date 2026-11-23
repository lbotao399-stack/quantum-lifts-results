"""Classical boundary minors in the three-scalar source of the sixth glueball."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import cache
from itertools import combinations, product, permutations
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'global_boundaries'))
from global_source_boundaries import colour, minor
from global_source_scan import square_determinant

# Five columns of length three in the ordinary polarization of (3,3,3,3,3).
DEGREE=5
FULL=(1<<(3*DEGREE))-1
PAIRS=tuple(combinations(range(3),2))
CUTS=((3,3,0,0,0),(3,3,0,0,0),(3,3,0,0,0),(3,3,1,0,0),
      (3,3,1,0,0),(3,2,1,0,0),(3,3,2,0,0),(3,2,2,0,0))

def parity(x):
    return (-1)**sum(a>b for j,a in enumerate(x) for b in x[j+1:])

@cache
def subsets(cuts):
    choices=[list(combinations(range(3),k)) for k in cuts]
    return [(sum(1<<(3*g+i) for g,s in enumerate(ch) for i in s),
             (-1)**sum(i-j for t in ch for j,i in enumerate(t))) for ch in product(*choices)]

@cache
def extensions(token,mask):
    role,*cols=token
    values=[[j for j in range(3) if mask & (1<<(3*c+j))] for c in cols]
    result=[]
    for labels in product(*values):
        if len(set(labels)) != len(labels):
            continue
        prev=mask
        sign=parity(labels)
        for c,i in reversed(tuple(zip(cols,labels))):
            prev ^= 1<<(3*c+i)
            sign*=(-1)**(((prev>>(3*c+i+1))&((1<<(2-i))-1)).bit_count())
        key=labels[0] if role=='x' else 9 if role=='j' else (3 if role=='y' else 6)+PAIRS.index(tuple(sorted(labels)))
        result.append((prev,sign,key))
    return result

@cache
def derivative(word):
    result=[]; phase=1
    for i,token in enumerate(word):
        if isinstance(token,str):continue
        tag,*c=token
        if tag=='y':
            result.append((phase,word[:i]+(('k',*c),)+word[i+1:]))
            phase=-phase
        elif tag=='z':
            result.append((phase,word[:i]+(('j',*c),)+word[i+1:]))
    return result

def bank(rng,cuts,marks,count):
    result=[]
    for index in range(count):
        cols=list(map(int,rng.permutation([j for j,k in enumerate(cuts) for _ in range(k)])))
        sizes=[1]*len(cols) if marks==0 else [2]+[1]*(len(cols)-2) if marks==1 else (
            [3]+[1]*(len(cols)-3) if index%2==0 else [2,2]+[1]*(len(cols)-4))
        start=0;word=[]
        for n in rng.permutation(sizes):
            n=int(n);word.append((('x','y','z')[n-1],*cols[start:start+n]));start+=n
        if index%3:
            pos=int(rng.integers(2,len(word)))
            deco=('ad','jo','trace','transpose')[int(rng.integers(4))]
            word.insert(pos,deco)
        if index%5==0 and len(word)>6:
            word.insert(len(word)-2,'ad' if index%2 else 'jo')
        result.append(tuple(word))
    return result

@cache
def recipe(order,case):
    cut=CUTS[case];other=tuple(3-z for z in cut)
    rng=np.random.default_rng(579103+case+97*order)
    left_marks,right_marks=((1,0) if order==1 else (1,1) if case%3!=1 else (2,0))
    return cut,bank(rng,cut,left_marks,30),bank(rng,other,right_marks,36),bool(right_marks)

class Point:
    def __init__(self,i,p,closed=False):
        gens,self.frames=colour(p)
        self.p=p
        rng=np.random.default_rng(741411+i)
        jets=list(((rng.integers(-3,4,size=(6,52))@gens.reshape(52,676)).reshape(6,26,26)%p).astype(float))
        def comm(a,b):return (a@b-b@a)%p
        splits=[comm(jets[j],jets[k]) for j,k in PAIRS]
        if closed:jets[3:6]=splits
        jac=(comm(jets[0],jets[5])-comm(jets[1],jets[4])+comm(jets[2],jets[3]))%p
        self.jets=jets+splits+[jac]
        self.eye=np.eye(26)
        self.zero=np.zeros_like(self.eye)
        from marked_jet_blocks import invariant_metric
        self.metric,self.inv=invariant_metric(p)
    def values(self):
        @cache
        def calculate(word,mask):
            if not word:
                assert mask==0
                return self.eye
            token=word[-1];prefix=word[:-1]
            if isinstance(token,str):
                x=calculate(prefix,mask)
                if token=='trace':return self.eye*(x.trace()%self.p)
                if token=='transpose':return (self.inv@x.T%self.p)@self.metric%self.p
                primal,dual,gram=self.frames[token]
                return ((x.reshape(676)@dual%self.p@gram%self.p)@primal%self.p).reshape(26,26)
            terms=extensions(token,mask)
            if not terms:return self.zero
            left=np.concatenate([calculate(prefix,m) for m,_,_ in terms],axis=1)
            right=np.concatenate([self.jets[j]*s for m,s,j in terms])
            return left@right%self.p
        return calculate
    def evaluate(self,order,case):
        cut,a,b,two=recipe(order,case);read=self.values();p=self.p
        w=len(a)*len(b);out=np.zeros((4,len(a),len(b)))
        def rows(words,mask,diff):
            if diff:return np.stack([sum((s*read(v,mask) for s,v in derivative(t)),self.zero)%p for t in words])
            return np.stack([read(t,mask) for t in words])
        def pair(a,b):
            x=a.reshape(len(a),676);y=b.reshape(len(b),676)
            vals=[x@b.transpose(0,2,1).reshape(len(b),676).T,
                  np.trace(a,axis1=1,axis2=2)[:,None] * np.trace(b,axis1=1,axis2=2)]
            for channel in ('ad','jo'):
                _,dual,gram=self.frames[channel]
                vals.append((x@dual%p@gram%p)@(y@dual%p).T)
            return np.asarray(vals)%p
        for mask,sg in subsets(cut):
            n=FULL^mask
            raw=pair(rows(a,mask,True),rows(b,n,False))
            if two:raw-=pair(rows(a,mask,False),rows(b,n,True))
            out+=sg*raw
            read.cache_clear()
        return np.asarray(out.reshape(4*w)%p,dtype=np.uint16)

@cache
def chosen(order,case):
    with np.errstate(all='ignore'):
        row=Point(35911,1009).evaluate(order,case)
    live=np.flatnonzero(row)
    rng=np.random.default_rng(432741+order*200+case)
    return np.sort(rng.choice(live,min(480,len(live)),replace=False))

def point(args):
    i,p,order,cases=args
    with np.errstate(all='ignore'):
        e=Point(i,p)
        return np.concatenate([e.evaluate(order,c)[chosen(order,c)] for c in cases])

def run(order,p,workers,points,dest):
    cases=list(range(len(CUTS)))
    for c in cases:chosen(order,c)
    limit=(429 if order==1 else 1863)
    assert 1 <= points <= limit
    ids=range(points)
    raw=np.empty((points,sum(len(chosen(order,c)) for c in cases)),dtype=np.uint16)
    checkpoints={min(points,n) for n in (192,420,840,1440,points)}
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        for k,v in enumerate(pool.map(point,((i,p,order,cases) for i in ids),chunksize=2),1):
            raw[k-1]=v
            if k in checkpoints:
                with np.errstate(all='ignore'):rows,cols,det=minor(raw[:k],p)
                record=dict(order=order,prime=p,points=k,rows=rows,columns=cols,
                            determinant=det,rank=len(rows),widths=[len(chosen(order,c)) for c in cases])
                output=dest/f'three-{order}-{p}-{k}.json'
                tmp=output.with_suffix('.part');tmp.write_text(json.dumps(record)+'\n');tmp.replace(output)
                print(order,p,k,len(rows),det,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['search','replay','controls'])
    ap.add_argument('--order',type=int,choices=[1,2],default=2)
    ap.add_argument('--prime',type=int,default=1009,choices=[1009,1013])
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--points',type=int,default=16)
    ap.add_argument('--out',type=Path,default=Path.cwd())
    ap.add_argument('--certificate',type=Path)
    args=ap.parse_args()
    if args.action=='search':
        run(args.order,args.prime,args.workers,args.points,args.out)
    elif args.action=='replay':
        cert=json.loads(args.certificate.read_bytes())
        assert cert['prime']==1009 and cert['order']==args.order and args.prime==1013
        cases=list(range(len(CUTS)))
        for c in cases: chosen(args.order,c)
        with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
            jobs=((i,args.prime,args.order,cases) for i in cert['rows'])
            raw=np.asarray([v[cert['columns']] for v in pool.map(point,jobs,chunksize=2)])
        with np.errstate(all='ignore'):
            det=square_determinant(raw,args.prime)
        result=dict(order=args.order,prime=args.prime,size=len(raw),determinant=det)
        (args.out/f'three-replay-{args.order}.json').write_text(json.dumps(result)+'\n')
        print(result,flush=True)
    else:
        with np.errstate(all='ignore'):
            e=Point(7891,args.prime,True)
            for c in range(len(CUTS)):
                assert not e.evaluate(2,c).any()
            positive=Point(7891,args.prime).evaluate(1,0)
            assert positive.any()
        print('All five-column Jacobi and Koszul boundary controls passed.')
