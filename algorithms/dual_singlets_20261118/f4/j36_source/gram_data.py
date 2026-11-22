"""Exact inputs for the last nine-dimensional Gram block in the F4 endpoint."""
from fractions import Fraction
from functools import cache
import io
from itertools import combinations, permutations
import json
from math import prod
from pathlib import Path
import tarfile

from endpoint_metric import ScalarLift, read, sign
from native import matrix
from verify import contents


def dependency(z, filename, member):
    blob = contents(z, 'independent-audit/dependencies/'+filename+'.tar.gz')
    with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as t:
        return t.extractfile(member).read()


def export(z, target):
    colour=read(z,'collected-O0/color-basis-native.json')
    lift=ScalarLift(colour)
    metric=matrix(colour['metric'])
    rows=json.loads(dependency(z,'f4-j36-S6-E2-exact','f4_j36_s6/f4-a5-prv-K.json'))
    frame=next(r['basis'] for r in rows if r['i']==r['j']==1)
    longest=json.loads(dependency(z,'f4-j36-native-active-projected-q-exact16-and-four-closed-cochains',
                           'f4_j36_active_project_exact/f4-j30-w0-native.json'))
    @cache
    def wedge(a,b):
        return sum(sign(p)*prod(metric.get((a[i],b[j]),0) for i,j in enumerate(p))
                   for p in permutations(range(len(a))))
    highest=(19,29,31,45,46)
    lowest=tuple(next(j for i,j in metric if i==k) for k in highest)
    coordinates=tuple((tuple(t[:3]),t[3],i) for i,block in enumerate(frame) for t in block)
    triplets=[t for t in combinations(range(52),3)
              if all(sum(lift.weight[a][j] for a in t)==0 for j in range(4))]
    seed={}
    for t in triplets:
        v=[0]*9
        for u,c,i in coordinates:
            v[i]+=c*wedge(u,t)
        seed[t]=[sign(lowest)*x for x in v]
    # Twice the Weyl-rotated inverse metric has integer coefficients.
    dual=[[] for _ in range(52)]
    for (i,j),v in lift.inverse.items():
        new,s=longest[i]
        dual[j].append((new,2*Fraction(s)*v))
    r=tuple(longest[i][0] for i in sorted(lowest))
    assert tuple(sorted(r))==highest
    gamma=sign(r)*prod(Fraction(longest[i][1]) for i in lowest)
    scale=629356*gamma/wedge(highest,tuple(sorted(lowest)))
    assert scale==-157339 and len(triplets)==236
    roots=[]
    for i,w in enumerate(lift.weight):
        if w>(0,0,0,0):
            j=next(k for k in range(4,52) if lift.weight[k]==tuple(-x for x in w))
            roots.append((i,j,2/metric[i,j]))
    def write(*fields):
        assert all(Fraction(x).denominator==1 for x in fields)
        print(*(int(x) for x in fields),file=target)
    for w in lift.weight:
        write(*w)
    for i in range(4):
        write(*(2*lift.inverse.get((i,j),0) for j in range(4)))
    write(*highest)
    write(*sorted(lowest))
    write(scale)
    write(len(roots))
    for root in roots:
        write(*root)
    write(len(lift.bracket))
    for (i,j),terms in sorted(lift.bracket.items()):
        write(i,j,len(terms),*(z for k,v in terms.items() for z in (k,v)))
    for d in dual:
        write(len(d),*(z for k,v in d for z in (k,v)))
    write(len(seed))
    for t,v in seed.items():
        write(*t,*v)
