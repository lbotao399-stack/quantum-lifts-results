"""Relative SU(2) BMN tree cohomology by restriction to an invertible scalar frame."""
from argparse import ArgumentParser
from collections import defaultdict
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import gcd
from pathlib import Path
import json

W=(1,0,-1)  # E,H,F in the defining sl(2) basis
ZERO=(0,0,0)

def add_triple(s,t): return tuple(a+b for a,b in zip(s,t))

def unit(i):return tuple(int(k==i) for k in range(3))

def wedge(s,t):
    if s&t:return 0
    inversions=0
    while t:
        low=t&-t;t-=low;inversions+=(s//(2*low)).bit_count()
    return (-1)**inversions


def mul(left,right):
    out=defaultdict(int)
    for (m,u),a in left.items():
        for (n,v),b in right.items():
            sign=wedge(m,n)
            if sign:out[m|n,add_triple(u,v)]+=a*b*sign
    return {k:v for k,v in out.items() if v}

ONE={(0,ZERO):1}
PSI=tuple(tuple({(1<<(3*i+j),ZERO):1} for j in range(3)) for i in range(3))
AUX=tuple({(0,unit(i)):1} for i in range(3))

def times(factors):
    value=ONE
    for factor in factors:
        value=mul(value,factor)
    return value


def determinant(columns):
    result=defaultdict(int)
    for perm in permutations(range(3)):
        sg=(-1)**sum(a>b for i,a in enumerate(perm) for b in perm[i+1:]);t=times(col[row] for col,row in zip(columns,perm))
        for k,v in t.items(): result[k]+=sg*v
    return {k:v for k,v in result.items() if v}

COL=tuple(tuple(PSI[i][j] for i in range(3)) for j in range(3))
DETS=tuple(determinant([COL[i] for i in ids]) for ids in combinations_with_replacement(range(3),3))
AUX_DETS=tuple(determinant([AUX,COL[i],COL[j]]) for i,j in combinations_with_replacement(range(3),2))
PAIRS=tuple({m:v for i in range(3) for m,v in mul(AUX[i],COL[j][i]).items()} for j in range(3))


def powers(n):
    for i in range(n+1):
        for j in range(n-i+1):yield (i,j,n-i-j)


def family(num_phi,k,kind):
    direction,num_aux=kind
    if direction=='upper':
        surplus=(num_phi+num_aux-k)//3
        assert 3*surplus==num_phi+num_aux-k
        pairs=max(0,num_aux-max(0,surplus));det=max(0,-surplus)
        for pair_idx in combinations(range(3),pairs):
            for det_idx in combinations(range(len(DETS)),det):
                yield times([PAIRS[r] for r in pair_idx]+[DETS[r] for r in det_idx]),k-pairs-3*det,num_aux-pairs
    else:
        surplus=(num_phi-k-num_aux)//3
        assert 3*surplus==num_phi-k-num_aux
        det=max(0,-surplus)
        for with_aux in range(1+min(det,num_aux)):
            for det_idx in combinations(range(len(DETS)),det-with_aux):
                core=times([DETS[r] for r in det_idx]);multi=combinations_with_replacement(range(len(AUX_DETS)),with_aux)
                for ids in multi:
                    tail=times(AUX_DETS[r] for r in ids)
                    yield mul(core,tail),k-3*det+with_aux,num_aux-with_aux


@cache
def maskweight(mask):
    return sum(W[i]+W[j] for i in range(3) for j in range(3) if mask&(1<<(3*i+j)))

def weight(t):return sum(a*b for a,b in zip(t,W))

class Rank:
    def __init__(self,p):self.p=p;self.rows={}
    def insert(self,vector):
        p=self.p
        if p==0:
            v={k:b for k,b in vector.items() if b}
            while v:
                first=min(v);c=v[first]
                if first not in self.rows:
                    content=gcd(*v.values())*(1 if c>0 else -1)
                    self.rows[first]={k:b//content for k,b in v.items()};return True
                pivot=self.rows[first];g=gcd(c,pivot[first]);a=pivot[first]//g; b=c//g
                v={k:n*a for k,n in v.items()}
                for k,n in pivot.items():
                    new=v.get(k,0)-b*n
                    if new:v[k]=new
                    else:v.pop(k,None)
            return False
        v={k:b%p for k,b in vector.items() if b%p}
        while v:
            first=min(v);c=v[first]
            if first not in self.rows:
                inv=pow(c,-1,p);v={k:(b*inv)%p for k,b in v.items()};self.rows[first]=v;return True
            for k,b in self.rows[first].items():
                n=(v.get(k,0)-c*b)%p
                if n:v[k]=n
                elif k in v:del v[k]
        return False
    def __len__(self):return len(self.rows)


def coefficients(num_phi,k,kind,p):
    by_weight=defaultdict(lambda:Rank(p));usign=-1 if kind[0]=='upper' else 1
    for core,left,naux in family(num_phi,k,kind):
        if not core or left<0:continue
        for subset in combinations(range(9),left):
            mask=sum(1<<b for b in subset)
            ext=mul(core,{(mask,ZERO):1})
            for u in powers(naux):
                row={(n,add_triple(v,u)):b for (n,v),b in ext.items()}
                if row:
                    values={maskweight(n)+usign*weight(v) for n,v in row};assert len(values)==1
                    by_weight[values.pop()].insert(row)
    return {w:tuple(r.rows.values()) for w,r in by_weight.items()}

# Standard component conventions Q psi_1=[phi^2,phi^3], cyclically, Q f=[phi^I,psi_I].
BRACKET={(0,1):(0,-2),(0,2):(1,1),(1,2):(2,-2)}
BRACKET.update({(j,i):(mu,-v) for (i,j),(mu,v) in list(BRACKET.items())})
Q_PSI={}
for i in range(3):
    mu,v=BRACKET[(i+1)%3,(i+2)%3];Q_PSI[3*i+mu]=v
Q_F=tuple({3*i+n:v for (i,n),(mu,v) in BRACKET.items() if mu==j} for j in range(3))
assert all(sum(c*Q_PSI.get(i,0) for i,c in row.items())==0 for row in Q_F)

@cache
def on_monomial(mask,f):
    result={}
    for i,coefficient in Q_PSI.items():
        bit=1<<i
        if mask&bit: result[mask^bit,f]=(-1)**((mask&(bit-1)).bit_count())*coefficient
    # Place Qf on the left of the exterior monomial, after the even f's.
    for j,power in enumerate(f):
        if not power:continue
        n=tuple(v-(h==j) for h,v in enumerate(f))
        for i,c in Q_F[j].items():
            bit=1<<i;sign=wedge(bit,mask)
            if sign:result[mask|bit,n]=sign*power*c
    return result


def image(poly,f):
    result=defaultdict(int)
    for (mask,u),coef in poly.items():
        for (n,g),b in on_monomial(mask,f).items(): result[n,u,g]+=coef*b
    return {k:v for k,v in result.items() if v}


def source_kind(w):
    if w[1]==w[2]:return 'lower',w[0]-w[1]
    if w[0]==w[1]:return 'upper',w[1]-w[2]
    raise ValueError(w)


def run(w,p,ylo=None,yhi=None):
    d=sum(w);kind=source_kind(w);rows={}
    for y in range(1,d+1):
        if ylo is not None and y<ylo or yhi is not None and y>yhi:continue
        dims=[0,0];ranks=[Rank(p),Rank(p)]
        for fdegree in range(d//3+1):
            k=d-y-2*fdegree;phi=y-k-fdegree
            if not (0<=k<=9 and phi>=0):continue
            coeff=coefficients(phi,k,kind,p)
            for f in powers(fdegree):
                for charge in range(2):
                    for vec in coeff.get(charge-weight(f),()):
                        dims[charge]+=1;ranks[charge].insert(image(vec,f))
        if any(dims):rows[y]=dict(weight_dimensions=dims,weight_ranks=[len(r) for r in ranks],dimension=dims[0]-dims[1],out=len(ranks[0])-len(ranks[1]))
    for y,r in rows.items():
        if ylo is None or y>ylo:
            incoming=rows.get(y-1,{}).get('out',0); h=r['dimension']-r['out']-incoming;assert h>=0
            r.update(incoming=incoming,cohomology_upper=h)
    return rows

if __name__=='__main__':
    a=ArgumentParser();a.add_argument('weight',nargs=3,type=int);a.add_argument('--prime',type=int,default=1000003);a.add_argument('--start',type=int);a.add_argument('--end',type=int);a.add_argument('--output',required=True,type=Path);a=a.parse_args()
    if not __debug__:raise ValueError('assertions required')
    r=run(tuple(a.weight),a.prime,a.start,a.end)
    from su2_bmn_character50 import characters
    expected=characters([a.weight])[','.join(map(str,a.weight))];assert all(v['dimension']==expected.get(str(y),0) for y,v in r.items());assert all(int(k) in r for k in expected if (a.start is None or int(k)>=a.start)and(a.end is None or int(k)<=a.end))
    data=dict(weight=a.weight,prime=a.prime,by_length=r);print([(y,d.get('cohomology_upper')) for y,d in r.items() if d.get('cohomology_upper')]);a.output.write_text(json.dumps(data,indent=2)+'\n')
