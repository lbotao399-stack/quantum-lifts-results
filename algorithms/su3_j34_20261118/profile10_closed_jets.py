"""Test all three-fermion coefficients on exact classical zeroes of the projected SU(3) boundary.

The form is first evaluated on all 24 matter coordinates, then pulled back to the
rational kernel derived from Q_letter for this scalar jet. No cochain-span fit.
"""
import argparse, gzip, hashlib, json, random, sys, os
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np

HERE=Path(__file__).resolve().parent
RUNTIME=HERE.parent/'su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
sys.path.insert(0,str(RUNTIME))
import sp_environment
from fortuity.full.letters import Q_letter, letters_up_to
PARENT='9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14'
TARGET='cb4f49285a43a3390aab5ac1768bd362302c368834a2d280b9432f794f7ab9ca'

def code(l):return (2*l[0]+l[1])*8+sum(1<<i for i in l[2])
def kind(l):return len(l[2]),l[0]+l[1]
def odd(x):return (x&7).bit_count()==2
MV=(3,5,6) # sorted theta pairs (0,1), (0,2), (1,2)
Z=np.zeros((3,3),dtype=np.int64); I=np.eye(3,dtype=np.int64)
U=[]
for i in range(3):
 for j in range(3):
  if i!=j:
   a=Z.copy();a[i,j]=1;U.append(a)
U.extend([np.diag([1,0,-1]),np.diag([0,1,-1])]);U=np.array(U)

def rhs(source,projected):
    raw=gzip.decompress(source.read_bytes());assert hashlib.sha256(raw).hexdigest()==(TARGET if projected else PARENT)
    lines=raw.splitlines(keepends=True);assert len(lines)==(27780 if projected else 517468);terms=[];hh=hashlib.sha256()
    for line in lines:
        c,w=json.loads(line);w=tuple(tuple(map(code,t)) for t in w)
        if not projected and Counter(((x&7).bit_count(),(x//8).bit_count()) for t in w for x in t)!={(1,0):6,(1,1):2,(2,0):3}:continue
        hh.update(line);q=2*Fraction(c);assert q.denominator==1;terms.append((int(q),w))
    assert hh.hexdigest()==TARGET and len(terms)==27780
    return terms

def replacements(plan):
    p=json.loads(plan.read_bytes())['incoming'];assert p['mu']==[1,1,7,4,3] and p['target_length']==11
    pat=lambda s:Counter({tuple(k):n for k,n in s['profile']}); target=pat(p['targets'][10]);out={}
    for s in p['sources']:
        sc=pat(s)
        for l in letters_up_to(*p['mu'][:2]):
            if sc[kind(l)]<=0:continue
            terms=[]
            for c,x,y in Q_letter(l):
                h=sc.copy();h[kind(l)]-=1;h[kind(x)]+=1;h[kind(y)]+=1
                if +h==target:terms.append((c,code(x),code(y)))
            if terms:out[s['index'],code(l)]=terms
    assert sorted(Counter(k[0] for k in out).items())==[(7,1),(9,2),(16,3),(17,6)]
    assert sum(map(len,out.values()))==48
    return out

def rational_kernel(matrix):
    # Unique reduced row echelon form over Q. Free columns are in ascending order.
    a=[list(map(Fraction,row)) for row in matrix];m=len(a);piv=[]
    for j in range(24):
        i=len(piv);h=next((k for k in range(i,m) if a[k][j]),None)
        if h is None:continue
        a[i],a[h]=a[h],a[i];div=a[i][j];a[i]=[x/div for x in a[i]]
        for k in range(m):
            if k!=i and a[k][j]:
                t=a[k][j];a[k]=[x-t*y for x,y in zip(a[k],a[i])]
        piv.append(j)
    bas=[]
    for j in range(24):
        if j in piv:continue
        v=[Fraction(0)]*24;v[j]=Fraction(1)
        for i,k in enumerate(piv):v[k]=-a[i][j]
        n=lcm(*(x.denominator for x in v));ints=[int(n*x) for x in v];g=gcd(*ints)
        row=[x//g for x in ints];assert all(sum(x*y for x,y in zip(t,row))==0 for t in matrix);bas.append(row)
    assert len(piv)+len(bas)==24
    return np.array(bas,dtype=object).reshape(-1,24)

SCALARS=['repeated','nilpotent','mixednil','regular','cartan','jordan','principal','dualnil']
@lru_cache(None)
def scalar_tangents(tag):
    H=np.diag([1,1,-2]);N=U[0];D=np.diag([1,0,-1]);E13=U[1];E23=U[3]
    S={'repeated':(H,-H,H),'nilpotent':(N,-N,N),'mixednil':(N,E13,N-E13),
       'dualnil':(E13,E23,E13-E23),'regular':(D,-D,D),'cartan':(D,H,D+H),
       'jordan':(H,N,H+N),'principal':(N+E23,E13,N+E23+E13)}[tag]
    assert all(not np.any(x@y-y@x) for x,y in combinations(S,2))
    rows=[]
    for i,j in combinations(range(3),2):
        a=np.zeros((3,8,3,3),dtype=np.int64)
        a[i]=U@S[j]-S[j]@U;a[j]=S[i]@U-U@S[i];rows.extend(a.reshape(24,9).T.tolist())
    return S,rational_kernel(rows).astype(np.int64)

def ring_jet(case):
    tag,seed=case;r=random.Random(5933781+seed)
    S,basis=scalar_tangents(tag);fields=dict(zip((1,2,4),S))
    for shift in (8,16):
        row=np.array([r.choice((-1,0,1)) for _ in basis])@basis
        matrices=(row.reshape(3,8)@U.reshape(8,9)).reshape(3,3,3)
        for x,m in zip((1,2,4),matrices):fields[shift+x]=m
    return fields

def kernel(fields,repls):
    linear=[];bosons=[]
    for (s,_),ts in repls.items():
        if s in (16,17):
            val=sum((c*(fields[x]@fields[y]) for c,x,y in ts),Z.copy());assert not val.any();bosons.append(val)
        else:
            val=np.zeros((3,8,3,3),dtype=np.int64)
            for c,x,y in ts:
                if x in MV:val[MV.index(x)]+=c*(U@fields[y])
                else:val[MV.index(y)]+=c*(fields[x]@U)
            linear.extend(val.reshape(24,9).T.tolist())
    assert len(bosons)==9 and len(linear)==27
    K=rational_kernel(linear)
    # Independent ordered-product evaluation for every odd tangent vector, on all 12 letter rules.
    for row in K:
        f=dict(fields);mat=row.reshape(3,8)@U.reshape(8,9).astype(object);f.update(zip(MV,mat.reshape(3,3,3)))
        for terms in repls.values():assert not sum((int(c)*(f[x]@f[y]) for c,x,y in terms),Z.copy()).any()
    return K

def evaluate(fields,terms):
    # Treat only the bosonic segments as numerical. Odd slots retain their order in each trace.
    @lru_cache(None)
    def segment(xs):
        m=I
        for x in xs:m=m@fields[x]
        return m
    @lru_cache(None)
    def trace(tr):
        pos=[i for i,x in enumerate(tr) if odd(x)]
        if not pos:return int(np.trace(segment(tr)))
        i=pos[0];wrap=tr[(pos[-1]+1):]+tr[:i];key=tuple(tr[j] for j in pos)
        if len(pos)==1:a=np.einsum('aij,ji->a',U,segment(wrap))
        else:
            j=pos[1];a=np.einsum('aij,bji->ab',U@segment(tr[i+1:j]),U@segment(wrap));assert len(pos)==2
        return key,a
    # Store coefficients with a fixed ordered slot sequence, before antisymmetrizing.
    by=defaultdict(lambda:np.zeros((8,8,8),dtype=np.int64));bnd=Counter(); nz=0;counts=Counter()
    for c,w in terms:
        factors=[]; key=[]
        for t in w:
            r=trace(t)
            if isinstance(r,int):c*=r
            else:key.extend(r[0]);factors.append(r[1])
        if not c:continue
        k=tuple(key);assert len(k)==3 and len(factors)>=2
        bound=abs(c)
        for a in factors:bound*=int(np.max(np.abs(a)))
        if not bound:continue
        bnd[k]+=bound;assert max(abs(c),bnd[k]) < 2**62
        tensor=factors[0]
        for a in factors[1:]:tensor=np.multiply.outer(tensor,a)
        by[k]+=c*tensor;nz+=1;counts[k]+=1
    # int64 bounds for all intermediate *trace* products, in addition to per-bucket accumulation.
    # max row-sum norm for a scalar; each word contains six scalars and precisely two scalar jets.
    norm=lambda m:max(sum(abs(int(x)) for x in row) for row in m)
    bound=3*max(1,*(norm(fields[x]) for x in (1,2,4)))**6*max(1,*(norm(v) for k,v in fields.items() if k>7))**2
    assert bound<2**60 and 6*sum(bnd.values())<2**63
    T=np.zeros((24,24,24),dtype=np.int64)
    for k,arr in by.items():T[tuple(slice(8*MV.index(x),8*(MV.index(x)+1)) for x in k)]=arr
    A=T+T.transpose(1,2,0)+T.transpose(2,0,1)-T.transpose(1,0,2)-T.transpose(2,1,0)-T.transpose(0,2,1)
    return A,nz

def task(case):
    f=ring_jet(case);K=kernel(f,REPL);A,nz=evaluate(f,TERMS);m=len(K)
    # Contract over arbitrary-precision integers before extracting the independent wedge coordinates.
    AO=A.astype(object);first=np.tensordot(K,AO,axes=(1,0));second=np.tensordot(K,first,axes=(1,1)).transpose(1,0,2);T=np.tensordot(second,K,axes=(2,1))
    ix=list(combinations(range(m),3));coords=[int(T[i,j,k]) for i,j,k in ix];raw=json.dumps(coords,separators=(',',':')).encode();nonzero=next(((i,int(c)) for i,c in zip(ix,coords) if c),None)
    if ix:
        i,j,k=nonzero[0] if nonzero else ix[0];assert int(T[i,j,k])==sum(int(K[i,a])*int(K[j,b])*int(K[k,c])*int(A[a,b,c]) for a,b,c in zip(*np.nonzero(A)))
    return dict(kind=case[0],seed=case[1],scalar_tangent_dimension=len(scalar_tangents(case[0])[1]),odd_kernel_dimension=m,cubic_coordinates=len(coords),nonzero=sum(bool(c) for c in coords),
                first_nonzero=nonzero,coordinate_sha256=hashlib.sha256(raw).hexdigest(),twice_rhs_unrestricted_nonzero=int(np.count_nonzero(A)//6),
                active_trace_terms=nz,noncommuting_even_pairs=sum(bool(np.any(x@y-y@x)) for x,y in combinations(f.values(),2)),
                fields={str(k):v.tolist() for k,v in f.items()},kernel_basis=K.tolist())

def initialize(source,plan,projected):
    global TERMS, REPL
    TERMS=rhs(source,projected);REPL=replacements(plan)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rhs',type=Path,default=HERE/'negative_Q1.jsonl.gz');p.add_argument('--projected',action='store_true');p.add_argument('--plan',type=Path,default=HERE/'classical_order1.json');p.add_argument('--workers',type=int,default=1);p.add_argument('--seeds',type=int,default=1);p.add_argument('--kinds',nargs='+',default=SCALARS,choices=SCALARS);p.add_argument('--output',type=Path);args=p.parse_args();assert __debug__ and args.seeds>0
    initialize(args.rhs,args.plan,args.projected);cases=[(k,i) for i in range(args.seeds) for k in args.kinds]
    if args.workers==1: rec=list(map(task,cases))
    else:
        import multiprocessing as mp
        with ProcessPoolExecutor(max_workers=min(args.workers,len(cases)),mp_context=mp.get_context('fork')) as pool:rec=list(pool.map(task,cases))
    data=dict(parent_uncompressed_sha256=PARENT,profile10_uncompressed_sha256=TARGET,projected_rhs_terms=len(TERMS),
              ordered_replacements={str(s):sum(len(ts) for (i,l),ts in REPL.items() if i==s) for s in (7,9,16,17)},
              coefficient_convention='2 times -Q1 O34; full cubic pullback in integral bases of the rational kernels',results=rec)
    text=json.dumps(data,indent=2)+'\n';print([(r['kind'],r['seed'],r['odd_kernel_dimension'],r['nonzero'],r['twice_rhs_unrestricted_nonzero']) for r in rec],flush=True)
    if args.output:args.output.write_text(text)
if __name__=='__main__':main()
