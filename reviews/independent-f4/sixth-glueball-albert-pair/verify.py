"""Rebuild two classical minors in the selected F4 sixth-glueball source."""
from collections import Counter
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path
import gzip
import json

import numpy as np
from exceptional_colour import projectors

HERE=Path(__file__).resolve().parent
PAIRS=tuple(combinations(range(4),2))
SEEDS={"separate":(0,3,6,14,16,21,67),"same":(5,26,36,55),"short":(1,)}
OFFSETS={"separate":996331,"same":821433,"short":114163}


def sign(w):
    return (-1)**sum(x>y for i,x in enumerate(w) for y in w[i+1:])


@cache
def flags():
    blocks=[tuple((p,sign(p)) for p in permutations(range(n))) for n in (4,4,3,3)]
    return [(a+b+c+d,sa*sb*sc*sd) for (a,sa),(b,sb),(c,sc),(d,sd) in product(*blocks)]


def bracket(i,j):
    if i==j:return 0,0
    return (1 if i<j else -1),4+PAIRS.index(tuple(sorted((i,j))))


def polynomial(kind,seed,distinguish_mark=False):
    order=np.random.Generator(np.random.PCG64(OFFSETS[kind]+seed)).permutation(14).tolist()
    out=Counter()
    for flag,coefficient in flags():
        w=tuple(flag[i] for i in order)
        if kind=="short":
            if len(set(w[:3]))!=3:continue
            s,m=bracket(w[5],w[6])
            traces=(tuple(sorted(w[3:5])),tuple(sorted(w[:3])))
            q1=tuple(sorted((m,*w[7:10])));q2=tuple(sorted(w[10:14]))
            coefficient*=s*sign(w[:3])
        else:
            traces=tuple(sorted((tuple(sorted(w[:2])),tuple(sorted(w[2:4])))))
            if kind=="separate":
                s,m=bracket(w[4],w[5]);t,n=bracket(w[9],w[10])
                q1=tuple(sorted((m,*w[6:9])));q2=tuple(sorted((n,*w[11:14])))
            else:
                s,m=bracket(w[8],w[9]);t,n=bracket(w[10],w[11])
                q1=tuple(sorted(w[4:8]));q2=tuple(sorted((m,n,*w[12:14])))
            coefficient*=s*t
        if coefficient:
            arms=(q1,q2) if distinguish_mark else tuple(sorted((q1,q2)))
            out[traces,arms]+=coefficient
    return {k:v for k,v in out.items() if v}


def weights(indices):
    return Counter(tuple(w.count(i) for i in range(4)) for w in indices)


def tensor(a,b):
    out=Counter()
    for u,n in a.items():
        for v,m in b.items():out[tuple(x+y for x,y in zip(u,v))]+=n*m
    return out


def symmetric_square(a):
    out=tensor(a,a)
    for u,n in a.items():out[tuple(2*x for x in u)]+=n
    assert all(v%2==0 for v in out.values())
    return Counter({k:v//2 for k,v in out.items()})


def multiplicity(ch):
    highest=(4,4,4,2)
    return sum(sign(w)*ch[tuple(highest[i]+w[i]-i for i in range(4))] for w in permutations(range(4)))


def characters():
    v=weights(combinations(range(4),1))
    wedge=weights(combinations(range(4),2));wedge3=weights(combinations(range(4),3))
    sym2,sym3,sym4=(weights(combinations_with_replacement(range(4),n)) for n in (2,3,4))
    traces=symmetric_square(sym2); arm=tensor(wedge,sym3)
    pure=tensor(traces,sym4)
    short=tensor(sym2,tensor(arm,sym4))
    types={
       "separate":(traces,tensor(arm,arm),symmetric_square(arm)),
       "same":(tensor(pure,sym2),tensor(wedge,wedge),symmetric_square(wedge)),
       "short":(short,tensor(wedge,v),wedge3),
       "nested":(tensor(pure,sym3),tensor(wedge,v),tensor(wedge,v)-wedge3)}
    return {k:[multiplicity(tensor(common,a)),multiplicity(tensor(common,b))] for k,(common,a,b) in types.items()}


def representation():
    data=json.loads(gzip.decompress((HERE/"representation26.json.gz").read_bytes()))
    basis=np.zeros((52,26,26),dtype=np.int64)
    for i,entries in enumerate(data["matrices"]):
        for a,b,c in entries:basis[i,a,b]=c
    return basis


def column(number,p,basis,projection,inverse,polys):
    coefficients=np.random.Generator(np.random.PCG64(144336+number)).integers(-3,4,size=(4,52))
    x=(coefficients@basis.reshape(52,-1)).reshape(4,26,26)%p
    letters=list(x)+[(x[a]@x[b]-x[b]@x[a])%p for a,b in PAIRS]
    @cache
    def trace(t):
        if len(t)==2:return int(np.sum(x[t[0]]*x[t[1]].T))%p
        return int(np.sum(letters[4+PAIRS.index(t[:2])]*x[t[2]].T))%p
    @cache
    def arm(q):
        mat=letters[q[0]]
        for k in q[1:]:mat=mat@letters[k]%p
        return projection@mat.reshape(-1)%p
    @cache
    def joined(q,r):return int(((arm(q)@inverse)%p)@arm(r))%p
    answers=[sum(c*trace(t)*trace(u)*joined(q,r) for ((t,u),(q,r)),c in row.items())%p for row in polys]
    # Every 26 by 26 matrix product has magnitude at most 26(p-1)^2
    # per entry. The outer polynomial sum uses unbounded Python integers.
    trace.cache_clear();arm.cache_clear();joined.cache_clear()
    return answers


def determinant(matrix,p):
    a=[list(row) for row in matrix]; n=len(a); result=1
    for j in range(n):
        i=next((i for i in range(j,n) if a[i][j]%p),None)
        if i is None:return 0
        if i!=j:a[i],a[j]=a[j],a[i];result=-result
        pivot=a[j][j]%p;result=result*pivot%p
        for i in range(j+1,n):
            m=a[i][j]*pow(pivot,-1,p)%p
            a[i]=[(b-m*c)%p for b,c in zip(a[i],a[j])]
    return result%p


def calculate():
    chars=characters()
    assert chars=={"separate":[9,7],"same":[5,4],"short":[5,1],"nested":[3,3]}
    recipes=[(k,s) for k,ss in SEEDS.items() for s in ss]
    rows=[polynomial(k,s) for k,s in recipes]
    basis=representation();data={}
    for p in (1009,1013):
        proj,inv=projectors(HERE,basis,p)[12]
        columns=[column(n,p,basis,proj,inv,rows) for n in range(12)]
        matrix=[list(row) for row in zip(*columns)]
        det=determinant(matrix,p)
        assert det=={1009:413,1013:948}[p]
        data[str(p)]=dict(matrix=matrix,determinant=det)
    return dict(shape=[4,4,4,2],character_multiplicities=chars,recipes=recipes,
                points=list(range(12)),minors=data)


if __name__=="__main__":
    print(json.dumps(calculate(),indent=2))
