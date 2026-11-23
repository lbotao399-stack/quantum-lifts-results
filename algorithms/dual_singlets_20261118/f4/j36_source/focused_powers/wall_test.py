"""Cartan restrictions of closed one-mark sources at the two F4 glueball levels."""
from fractions import Fraction
from itertools import permutations
import gzip
import json
from pathlib import Path
from math import prod

def determinant(a):
    n=len(a)
    return sum((-1)**sum(i>j for k,i in enumerate(p) for j in p[k+1:])
               *prod(a[k][j] for k,j in enumerate(p)) for p in permutations(range(n)))

def inverse(a):
    return [[Fraction((-1)**(i+j)*determinant([[x for k,x in enumerate(r) if k!=i]
                                               for l,r in enumerate(a) if l!=j]),determinant(a))
             for j in range(len(a))] for i in range(len(a))]

def dot(x,y):return sum(a*b for a,b in zip(x,y))
def multiply(m,x):return tuple(dot(t,x) for t in m)
path=Path(__file__).resolve().parent.parent/'first_page'
raw=json.loads(gzip.decompress((path/'representation26.json.gz').read_bytes()))
diag=[[0]*26 for _ in range(4)]
for j,row in enumerate(raw['matrices'][:4]):
    for k,l,v in row:
        assert k==l
        diag[j][k]=v
weights=tuple(zip(*diag))
metric=[[Fraction(dot(x,y),6) for y in diag] for x in diag]
dual=inverse(metric)
def bilinear(x,y):return dot(x,multiply(metric,y))
def normal(triple):
    cov=tuple((-1)**(3+j)*determinant([[v for k,v in enumerate(row) if k!=j] for row in triple]) for j in range(4))
    return multiply(dual,cov)
def values(triple,marks):
    n=normal(triple)
    return (bilinear(n,n)**2*sum(bilinear(x,y) for x,y in zip(triple,marks)),
            sum(dot(w,n)**4*sum(dot(w,x)*dot(w,y) for x,y in zip(triple,marks)) for w in weights))
zero=(0,)*4

def run():
    native=json.loads(gzip.decompress((path/'inputs.json.gz').read_bytes()))
    brackets={(i,j):dict(v) for i,j,v in native["brackets"]}
    rows=[]
    for label in native['raising']:
        root=tuple(brackets.get((i,label),{}).get(label,0) for i in range(4))
        if not any(root):raise ValueError('zero root')
        h=multiply(dual,root)
        pivot=next(i for i,x in enumerate(root) if x)
        plane=tuple(tuple((root[pivot] if j==i else 0)-(root[i] if j==pivot else 0) for j in range(4)) for i in range(4) if i!=pivot)
        assert len(plane)==3 and all(dot(root,x)==0 for x in plane)
        n=normal(plane)
        totals=[0,0]
        for i in range(3):
            changes=list(plane)
            changes[i]=h
            dn=normal(changes)
            # The mark in the positive cyclic pair complementary to i has colour h.
            a=plane[i]
            totals[0]+=4*bilinear(n,n)*bilinear(n,dn)*bilinear(a,h)+bilinear(n,n)**2*bilinear(h,h)
            totals[1]+=sum(4*dot(w,n)**3*dot(w,dn)*dot(w,a)*dot(w,h)
                           +dot(w,n)**4*dot(w,h)**2 for w in weights)
            # Independent coefficient extraction; changing one scalar makes degree at most five.
            samples=[]
            for t in range(6):
                test=list(plane);test[i]=tuple(u+t*v for u,v in zip(a,h))
                samples.append(values(test,tuple(h if j==i else zero for j in range(3))))
            weights_d=tuple(Fraction(n,d) for n,d in [(-137,60),(5,1),(-5,1),(10,3),(-5,4),(1,5)])
            coeff=[dot(weights_d,col) for col in zip(*samples)]
            if i==0:direct=[0,0]
            direct=[a+b for a,b in zip(direct,coeff)]
        assert direct==totals
        rows.append(dict(root=label,length=str(bilinear(h,h)),equations=list(map(str,totals))))
    return rows

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
