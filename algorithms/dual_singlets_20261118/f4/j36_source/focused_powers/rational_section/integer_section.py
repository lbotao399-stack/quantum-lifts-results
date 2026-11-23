"""Integral interpolation of the original 703-pivot section and the ordered crossed graph."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from math import gcd
from multiprocessing import get_context
from pathlib import Path
import gzip
import json
import os
import numpy as np
try:
    from scipy.sparse import csr_matrix
except ImportError:
    csr_matrix = None
import sys

FIRST_PAGE = Path(__file__).resolve().parents[2] / "first_page"
sys.path.insert(0, str(FIRST_PAGE))
from double_trace_flag import Evaluation, family
from triple_trace_flag import Triple
from single_trace_flag import catalogue, native26
from flag_boundaries import trace_words

HERE=Path(__file__).resolve().parent
RANGES={'four':(((0,1,2),(0,2,3),(1,2,0),(2,0,1),(3,0,2)),
                ((0,1,2),(0,1,3),(0,2,3),(1,2,3),(0,0,1),(0,1,0))),
        'nested':(((0,1,2),(0,2,3),(2,0,1),(1,0,3)),
                  ((0,1,2),(0,1,3),(1,2,3),(0,0,1)))}
TRAIN=json.loads((HERE/'selected-1009.json').read_bytes())
PAIRS=((0,1),(0,2),(1,3))
CROSSED=((1,0,3),(0,1,3),(0,2,3),(1,2,3))


def bound(a):
    return int(max(abs(a.min()), abs(a.max())))


@cache
def words():
    table=json.loads((FIRST_PAGE/'triple_selected.json').read_bytes())['rows']
    names=[table[j] for kind,j in TRAIN['rows'] if kind=='reference']
    ordinary=[s for s in names if not s[0].startswith('T')]
    triples=[s for s in names if s[0].startswith('T')]
    pos={tuple(v):i for i,v in enumerate(ordinary+triples)}
    return Evaluation(ordinary),Triple(triples),[pos[tuple(v)] for v in names]


def powers(x):
    x=x.astype(np.float64);dim=x.shape[-1];p={1:x}
    for n in range(2,8):
        last=p[n-1]
        assert dim*bound(x)*bound(last)<2**52
        value=np.empty((4**n,dim,dim))
        for j in range(4):
            value[j::4]=(last.reshape(-1,dim)@x[j]).reshape(-1,dim,dim)
        p[n]=value
    return p


def traces(table,codes,n):
    a,b=table[n//2],table[n-n//2]
    ceiling=a.shape[1]**2*bound(a)*bound(b)
    out=np.empty(len(codes),dtype=np.int64 if ceiling<2**62 else object)
    for j in range(0,len(codes),128):
        keys=codes[j:j+128]
        left,right=a[keys//len(b)],b[keys%len(b)]
        if ceiling<2**53:
            part=np.einsum('nij,nji->n',left,right).astype(np.int64)
        elif ceiling<2**62:
            part=np.einsum('nij,nji->n',left.astype(np.int64),right.astype(np.int64))
        else:
            # A dense control can exceed the coarse whole-table bound.
            part=np.zeros(len(keys),dtype=object)
            for h,l in enumerate(limbs(left)):
                for k,r in enumerate(limbs(right)):
                    small=np.einsum('nij,nji->n',l,r).astype(np.int64)
                    part+=small.astype(object)*4096**(h+k)
        out[j:j+128]=part
    return out


def sum_groups(co,vals,index=None,starts=None):
    v=vals if index is None else vals[index]
    # Object arithmetic only when the absolute upper bound requires it.
    if v.dtype==object or int(np.abs(co).sum())*bound(v)>=2**62:
        v=v.astype(object)
    terms=co*v
    return sum(map(int,terms)) if starts is None else np.add.reduceat(terms,starts)


def multiply(a,b):
    if a.dtype==object or b.dtype==object or bound(a)*bound(b)>=2**62:
        return a.astype(object)*b
    return a*b


def trace_row(coef):
    one,three,ix=words()
    data={};left={};right={};triples=None
    for name in ('26','52'):
        rep=family(name);dim=rep.shape[1]
        table=powers((coef@rep.reshape(52,-1)).reshape(4,dim,dim))
        if name in one.names:
            co,col,starts=one.single
            data[name]=sum_groups(co,traces(table,catalogue(),14),col,starts)
        if name=='52' and one.short_terms:
            _,sets,size=trace_words();frame=np.zeros(size,dtype=np.int64);k=1
            for n,ids in sets.items():
                frame[k:k+len(ids)]=traces(table,ids,n);k+=len(ids)
            ans=[]
            for co,idx in one.short_terms:
                v=np.ones(idx.shape[1],dtype=np.int64)
                for jj in idx:v=multiply(v,frame[jj])
                ans.append(sum_groups(co,v))
            data['short']=ans
        for n,dd in one.doubles.items():
            left[n,name]=traces(table,dd.long,n)
            right[n,name]=traces(table,dd.short,14-n)
        if name=='26':
            vv={n:traces(table,t,n) for n,t in three.words.items()}
            triples=[0]*three.count
            for selected,starts,parts,ids,co in three.blocks:
                if not len(selected):continue
                v=np.ones(len(parts[0][1]),dtype=np.int64)
                for n,j in parts:v=multiply(v,vv[n][j])
                for j,val in zip(selected,sum_groups(co,v,ids,starts)):triples[j]=int(val)
        del table
    for name,(n,a,b) in one.pairs.items():
        dd=one.doubles[n]
        v=multiply(left[n,a][dd.left],right[n,b][dd.right]);out=[0]*dd.size
        for j,val in zip(dd.active,sum_groups(dd.coefficients,v,dd.indices,dd.starts)):
            out[j]=int(val)
        data[name]=out
    all_rows=[int(data[name][k]) for (name,_),k in zip(one.rows,one.positions)]+triples
    return [all_rows[k] for k in ix]


@cache
def tensors():
    raw=json.loads(gzip.decompress((FIRST_PAGE/'representation26.json.gz').read_bytes()))
    h=np.zeros((26,26),dtype=np.int64)
    for i,j,value in raw['form']:
        q=Fraction(value);assert q.denominator==1;h[i,j]=q.numerator
    d=3*h
    d[:2,:2]=[[2,1],[1,2]]
    assert np.array_equal(h@d,3*np.eye(26,dtype=np.int64))
    obj=json.loads((HERE/'cubic.json').read_bytes());assert obj['k']=='28/3'
    c=np.zeros((26,)*3,dtype=np.int64)
    for i,j,k,v in obj['cubic']:c[i,j,k]=v
    def corner(x,den):
        x=x.reshape(26,-1)
        assert 26**2*bound(x)**2*bound(d)<2**62
        z=((x.T@d)@x).reshape((26,)*4)
        g=gcd(den,int(np.gcd.reduce(z.reshape(-1))))
        return z//g,den//g
    low,dl=corner(c,28)
    hi,dh=corner((d@c)@d,28*81)
    exc,de=corner(d@c,28*9)
    return low.astype(float),hi.transpose(0,2,1,3).astype(float),exc.transpose(0,2,1,3).reshape(676,676),(dl,dh,de)


def exact_apply(mat,t,axes):
    front=tuple(range(len(axes)));old=np.moveaxis(t,axes,front)
    width=26**len(axes)
    # Each CSR row has the stated 1-norm, bounding every partial sum.
    m,norm=mat
    assert norm*bound(t)<2**52
    with np.errstate(all='ignore'):
        v=m@old.reshape(width,-1)
    return np.moveaxis(v.reshape((26,)*4),front,axes)


def action(a):
    norm = int(np.abs(a).sum(axis=1).max())
    if csr_matrix is not None:
        return csr_matrix(a,dtype=np.float64),norm
    return a.astype(np.float64),norm


SIGNS={p:(-1)**sum(i>j for k,i in enumerate(p) for j in p[k+1:])
       for n in (3,4) for p in permutations(range(n))}


def halves(x,t,arms):
    matrices=[action(a) for a in x]
    # Transposed Lie actions use the transpose of the native commutator.
    # That sign is supplied by the caller for the left half.
    marks={(i,j):action(x[i]@x[j]-x[j]@x[i]) for i,j in combinations(range(4),2)}
    ends={}
    for a,b,c in sorted({item[0] for item in arms}):
        total=np.zeros_like(t)
        for (i,j),m in marks.items():
            k,l=(r for r in range(4) if r not in (i,j))
            u=exact_apply(m,t,(a,))
            v=exact_apply(matrices[l],exact_apply(matrices[k],u,(b,)),(c,))
            w=exact_apply(matrices[k],exact_apply(matrices[l],u,(b,)),(c,))
            assert bound(total)+bound(v)+bound(w)<2**53
            total+=SIGNS[(i,j,k,l)]*(v-w)
        ends[a,b,c]=total
    out=[]
    for end,mid in arms:
        t0=ends[end];total=np.zeros_like(t0)
        for per in permutations(range(3)):
            u=t0
            for i,axis in zip(per,mid):u=exact_apply(matrices[i],u,(axis,))
            assert bound(total)+bound(u)<2**53
            total+=SIGNS[per]*u
        out.append(total.ravel())
    return np.stack(out)


@cache
def plan():
    result=[]
    for tag,i in TRAIN['rows']:
        if tag=='reference':continue
        ends,mids=RANGES[tag]
        n=len(ends)*len(mids)*len(ends)*len(mids)
        a,b,c,d=list(product(ends,mids,ends,mids))[i%n]
        transfers=() if tag=='four' else (PAIRS[i//n],)
        result.append(((a,b),(c,d),transfers))
    result.append((CROSSED[:2],CROSSED[2:],((1,2),(0,1))))
    return result


def limbs(raw):
    a=raw.astype(np.int64)
    assert np.array_equal(raw,a)
    sign=np.sign(a);v=abs(a)
    while v.any():
        yield ((v%4096)*sign).astype(float)
        v=v//4096


def gram(left,right):
    # Split signed integers into base-4096 limbs for exact BLAS products.
    # Each dot sums 26^4 terms, comfortably below the binary64 integer limit.
    assert left.shape[1]*4095**2<2**53
    ld=list(limbs(left));rd=list(limbs(right))
    out=np.zeros((len(left),len(right)),dtype=object)
    for j,a in enumerate(ld):
        for k,b in enumerate(rd):
            raw=a@b.T
            out+=raw.astype(np.int64).astype(object)*4096**(j+k)
    return out


def graph_row(coef):
    matrices=native26()
    x=(coef@matrices.reshape(52,676)).reshape(4,26,26)
    lo,hi,ex,den=tensors();info=plan()
    ll=sorted({a for a,b,t in info});rr=sorted({b for a,b,t in info})
    a=-halves(x.transpose(0,2,1),lo,ll)
    b=halves(x,hi,rr)
    ie=action(ex);rows=[0]*len(info)
    scales=[den[0]*den[1]*den[2]**len(t) for _,_,t in info]
    for chain in sorted({t for _,_,t in info}):
        pos=[i for i,(l,r,t) in enumerate(info) if t==chain]
        ls=sorted({ll.index(info[i][0]) for i in pos})
        rs=sorted({rr.index(info[i][1]) for i in pos})
        bb=b[rs]
        for pair in chain:bb=np.stack([exact_apply(ie,v.reshape((26,)*4),pair).ravel() for v in bb])
        z=gram(a[ls],bb)
        for i in pos:rows[i]=int(z[ls.index(ll.index(info[i][0])),rs.index(rr.index(info[i][1]))])
    return rows,scales


def sample(i,dense=False):
    if dense:
        return np.random.default_rng(144336+TRAIN['points'][i]).integers(-3,4,size=(4,52))
    rng=np.random.default_rng(244116+i)
    x=np.zeros((4,52),dtype=np.int64)
    for row in x:
        ids=rng.choice(52,12,replace=False)
        row[ids]=rng.choice((-1,1),size=len(ids))
    return x


def work(args):
    i,mode=args
    coef=sample(i)
    with np.errstate(all='ignore'):
        t=trace_row(coef) if mode!='graph' else []
        g,d=graph_row(coef) if mode!='trace' else ([],[])
    return i,t+g


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('mode',choices=['trace','graph','all','control'])
    ap.add_argument('--start',type=int,default=0);ap.add_argument('--end',type=int,default=1)
    ap.add_argument('--workers',type=int,default=1);ap.add_argument('--out',type=Path,default=Path('.'))
    args=ap.parse_args();print('scales',tensors()[3],flush=True)
    if args.mode=='control':
        # Pinned dense evaluation can exceed the sparse no-overflow bound. Check at a fresh point
        # against the archived modular evaluator by supplying identical coefficients to it.
        for k in range(args.start,args.end):
            i,vals=work((k,'all'))
            print(i,len(vals),max(abs(x).bit_length() for x in vals),flush=True)
    else:
        if args.mode!='graph':words()
        with ProcessPoolExecutor(args.workers,mp_context=get_context('fork')) as pool:
            for i,row in pool.map(work,((i,args.mode) for i in range(args.start,args.end)),chunksize=1):
                dst=args.out/('%s-%04d.json'%(args.mode,i));tmp=dst.with_suffix('.tmp')
                tmp.write_text(json.dumps(row,separators=(',',':'))+'\n');tmp.replace(dst)
                if i%20==0:print('completed point',i,flush=True)
