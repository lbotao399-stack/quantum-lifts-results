"""First differential of a scalar-bilinear dressing of the degree-15 singleton Pfaffian."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction
from itertools import combinations,product
from math import factorial,prod
import json
from pathlib import Path
import sys

from mixed_ward import normalize


def add(x,y):return tuple(a+b for a,b in zip(x,y))
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def parity(x):return (x[0]+x[1])%2


def plane_log_kernel(n):
    from orthogonal_glueball_first_pair import PHI, single_quantum
    from mixed_incoming import so_canon
    res=[]
    for length in (3,5,7):
        orbits=defaultdict(list)
        for w in product(range(3),repeat=length):
            sign,key=so_canon((tuple(PHI[i] for i in w),))
            if sign:orbits[key[0]].append((sign,w))
        for word,orbit in orbits.items():
            coordinates=Counter()
            for factors,c in single_quantum(n,word).items():
                if any(sum(t)>4 for t in factors):continue
                assert 1<=len(factors)<=2
                tr=2**len(factors)*(-1)**((length-1)//2)
                if len(factors)==1:
                    coordinates[factors]+=tr*c
                else:
                    x,y=factors
                    coordinates[factors]+=tr*c*(2 if x==y else 1)
                    merged=add(x,y)
                    if max(merged[:2])<=1 and sum(merged)<=4:
                        coordinates[merged,]+=tr*c*(-1)**(x[1]*y[0])
            if coordinates:res.append((length,orbit,coordinates))
    return res


def matmul(x,y):
    a,b,c,d=x;e,f,g,h=y
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h


def make_class(core,u,v,log):
    q=Counter();b=Counter()
    e=[tuple(int(j==i+2) for j in range(5)) for i in range(3)]
    for i,j in product(range(3),repeat=2):
        ab=matmul(core[i],core[j])
        q[add(e[i],e[j])]+=Fraction(ab[0]+ab[3],2)
        b[add(e[i],e[j])]+=-2*u[i]*v[j]
    kernel=Counter()
    for length,orbit,coordinates in log:
        coef=0
        for s,w in orbit:
            a=(1,0,0,1)
            for i in w:a=matmul(a,core[i])
            coef+=s*(a[0]+a[3])
        for w,c in coordinates.items():kernel[w]+=c*Fraction(coef,2*length)
    q={w:c for w,c in q.items() if c};b={w:c for w,c in b.items() if c}
    kernel={w:c for w,c in kernel.items() if c}
    zero=(0,)*5
    @lru_cache(None)
    def reductions(site):
        options=[(site,1)]
        for t,c in q.items():
            rem=sub(site,t)
            if min(rem)>=0:options.append((rem,c))
        return options
    @lru_cache(None)
    def core_value(sites):
        fixed={i for i,t in enumerate(sites) if t!=zero and (t[:2]!=(0,0) or t not in q)}
        if len(fixed)>2:return 0
        total=0
        for size in (1,2):
            if size<len(fixed):continue
            pool=[i for i,t in enumerate(sites) if i not in fixed and t!=zero]
            for mix in combinations(pool,size-len(fixed)):
                live=tuple(sorted(fixed|set(mix)))
                coef=prod(q.get(t,1) for i,t in enumerate(sites) if i not in live)
                for opt in product(*(reductions(sites[i]) for i in live)):
                    residual=tuple(w for w,n in opt)
                    if zero in residual:continue
                    s,key=normalize(residual)
                    if s:total+=coef*prod(n for w,n in opt)*s*kernel.get(key,0)
        return total
    @lru_cache(None)
    def score(sites):
        total=0
        for i,t in enumerate(sites):
            for w,c in b.items():
                rem=sub(t,w)
                if min(rem)>=0:
                    # Pure scalar factors commute with every site.
                    s,key=normalize(sites[:i]+(rem,)+sites[i+1:])
                    if s:total+=s*c*core_value(key)
        return total
    return score


@lru_cache(None)
def undo_substitutions(w,last=5):
    if last<0:return ((w,1),)
    spin,flavour=divmod(last,3);prefix=0;out=Counter()
    for i,t in enumerate(w):
        if t[spin]:
            new=list(t);new[spin]=0;new[flavour+2]+=1
            s,key=normalize(w[:i]+(tuple(new),)+w[i+1:])
            if s:
                scale=(-1)**(prefix+sum(t[:spin]))*new[flavour+2]*s
                for v,c in undo_substitutions(key,last-1):out[v]+=scale*c
        prefix+=sum(t[:2])
    return tuple((v,c) for v,c in out.items() if c)


def lowering_rows(sites):
    """Image of both simple flavour lowerings in fixed labelled-plane coefficients."""
    rows=defaultdict(Counter)
    for j,w in enumerate(sites):
        for lower in (0,1):
            for i,t in enumerate(w):
                if not t[3+lower]:continue
                u=list(t);u[3+lower]-=1;u[2+lower]+=1
                s,highest=normalize(w[:i]+(tuple(u),)+w[i+1:])
                if s:rows[lower,highest][j]+=s*u[2+lower]
    return [[row.get(j,0) for j in range(len(sites))] for row in rows.values()]


def run(a):
    sys.path.insert(0,str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from flint import fmpq_mat
    from terminal_commutators import words
    from mixed_ward import kernel_rows
    data=json.loads(a.relations.read_text())
    basis=words((4,4,4,2,2))
    inds=data['source_plane_indices']
    sites=[basis[i] for i in inds]
    assert min(map(len,sites))==7 and len(sites)==346
    lower=lowering_rows(sites)
    fin=[r['planes'] for r in data['rows']]
    probe=kernel_rows(fmpq_mat(lower+fin).transpose())
    print('high planes',len(sites),'lowering and finite ranks',fmpq_mat(lower).rank(),fmpq_mat(lower+fin).rank(),flush=True)
    log=plane_log_kernel(13)
    print('connected scalar log words',len(log),flush=True)
    fields=sorted({v for t in sites for v,c in undo_substitutions(t)})
    print('BMN first target planes',len(fields),flush=True)
    records=[]
    for line in a.parameters.read_text().splitlines()[:a.count]:
        i,*terms=map(int,line.split());assert len(terms)==15
        polar=[(h,e,f,-h) for h,e,f in (terms[j:j+3] for j in (0,3,6))]
        u,v=terms[9:12],terms[12:15]
        score=make_class(polar,u,v,log)
        values={w:score(w) for w in fields}
        raw=[sum(c*values[v] for v,c in undo_substitutions(w)) for w in sites]
        projected=fmpq_mat([list(map(str,raw))])*probe.transpose()
        records.append(dict(seed=i,matrices=polar,u=u,v=v,raw=list(map(str,raw)),projected=list(map(str,projected.entries()))))
        print('source',i,'nonzero',sum(bool(x) for x in projected.entries()),'rank',fmpq_mat([r['projected'] for r in records]).rank(),flush=True)
    a.out.write_text(json.dumps(dict(planes=inds,lower_rank=fmpq_mat(lower).rank(),finite_plus_lower_rank=fmpq_mat(lower+fin).rank(),
                                     fields=len(fields),connected=len(log),rows=records),separators=(',',':'))+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--relations',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--parameters',type=Path,required=True)
    p.add_argument('--count',type=int,default=4)
    run(p.parse_args())
