"""Exact SU(2) replay. Builds the nonabelian model explicitly; uses shared Wick/jet routines.

Input: the single hash-pinned witness only. No saved Q images are used to
construct either side. Stored products are compared *after* reconstruction.
Run via tools/project.py replay su2-boundary; direct execution also works.
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
from functools import lru_cache
import hashlib
import io
import itertools
import json
import pickle
import sys
import zipfile

if not __debug__:
    raise SystemExit("Checks require enabled assertions (remove -O).")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'sources/engines'))
sys.path.insert(0, str(ROOT / 'sources/all-loop/product_kernels/code'))
import qengine
import product_kernels as PK


class PureSU2(qengine.Theory):
    """c=c^A t_A, b=b_A t^A. [H,E]=2E, [H,F]=-2F, [E,F]=H."""
    def __init__(self):
        self.P, self.R, self.degW = 0, qengine.Ring(0), 3
        colors = ('e12', 'e21', 'h1')
        self.odd = {'c'+s for s in colors}
        self.partners = {p+s:q+s for p,q in [('b','c'),('c','b')] for s in colors}
        self.base = {s:{} for s in self.partners}
        vertices = {}
        for i,j,k,v in [('h1','e12','e12',2),('h1','e21','e21',-2),('e12','e21','h1',1)]:
            for a,b,c,f in [(i,j,k,v),(j,i,k,-v)]:
                term = self.mul(self.v('c'+a), self.v('c'+b))
                self.base['c'+c] = self.add(self.base['c'+c], self.scale(term, Fraction(f,2)))
                term = self.mul(self.v('c'+a), self.v('b'+c))
                self.base['b'+b] = self.add(self.base['b'+b], self.scale(term, -f))
                word, sign = self.canon((('b'+c,0,0),('c'+a,0,0),('c'+b,0,0)))
                key = tuple(l[0] for l in word)
                vertices[key] = vertices.get(key,0) + sign*Fraction(f,2)
        self.vertices = list(vertices.items())

    def rotation(self, p, axis):
        out=defaultdict(int)
        for w, coeff in p.items():
            for n,(field,i,j) in enumerate(w):
                m=(i,j)[axis]
                if m:
                    di,dj=((-1,1),(1,-1))[axis]
                    key, sign=self.canon(w[:n]+((field,i+di,j+dj),)+w[n+1:])
                    out[key] += sign*m*coeff
        return self.fin(out)


# Portable version of the exact forest contraction in sources/engines/exact_star.py.
# Both named input vertices and all surviving jet positions retain source order.
@lru_cache(None)
def topology(n, ti):
    forest=PK.forests(n)[ti]
    edges, roots=PK.edges_of(forest)
    return tuple(edges[:-1]), tuple(roots), PK.aut_forest(forest)


@lru_cache(1<<15)
def product_contraction(T,n,ti,a,b,sa,sb):
    edges,roots,_=topology(n,ti); A=n+1; B=n+2
    slots=dict(zip((e for e in edges if e[1]==A),sa));slots.update(zip((e for e in edges if e[1]==B),sb))
    choices=[]
    for v in range(1,n+1):
        needed=[T.partners[(a if j==A else b)[slots[(i,j)]]] for i,j in slots if i==v]
        choices.append([(w,k) for w,k in T.vertices if all(w.count(s)>=needed.count(s) for s in set(needed))])
    out=defaultdict(int)
    for verts in itertools.product(*choices):
        word=tuple((i,g) for i,(w,_) in enumerate(verts,1) for g in w)+tuple((A,g,k) for k,g in enumerate(a))+tuple((B,g,k) for k,g in enumerate(b))
        p={word:1}
        for i,j in reversed(edges):
            p=T.wiedge(p,i,j,slots.get((i,j)))
            if not p:break
        pref=1
        for _,c in verts: pref*=c
        for w,c in p.items():out[w]+=pref*c
    return tuple(T.fin(out).items())


def restricted_derivative(word,positions,jet):
    p={word:1}
    for axis,order in enumerate(jet,1):
        for _ in range(order):
            q=defaultdict(int)
            for w,c in p.items():
                for i in positions:
                    u=list(w[i]);u[axis]+=1;q[w[:i]+(tuple(u),)+w[i+1:]]+=c
            p=q
    return p


def star(T,n,p,q):
    if n==0:return T.mul(p,q)
    out=defaultdict(int); A=n+1; B=n+2
    for wa,ca in p.items():
        for wb,cb in q.items():
            for ti in range(len(PK.forests(n))):
                edges,roots,aut=topology(n,ti); ext=[e for e in edges if e[1]>n]
                na=sum(j==A for _,j in ext); nb=len(ext)-na
                for sa in itertools.permutations(range(len(wa)),na):
                    for sb in itertools.permutations(range(len(wb)),nb):
                        cs=product_contraction(T,n,ti,tuple(u[0] for u in wa),tuple(u[0] for u in wb),sa,sb)
                        if not cs:continue
                        ia,ib=iter(sa),iter(sb); jets=tuple((wa[next(ia)] if j==A else wb[next(ib)])[1:] for _,j in ext)
                        for tags,c in cs:
                            for jj,k in PK.moments(n,ti,jets).items():
                                w=[];apos=[]
                                for t in tags:
                                    if len(t)==2:w.append((t[1],*jj[roots.index(t[0])]))
                                    elif t[0]==A:apos.append(len(w));w.append(wa[t[2]])
                                    else:w.append(wb[t[2]])
                                for term,m in restricted_derivative(tuple(w),apos,jj[-1]).items():
                                    key,sign=T.canon(term);out[key]+=sign*m*c*k*ca*cb*Fraction((-2)**n,aut)
    return T.fin(out)


class RationalWitness(pickle.Unpickler):
    def find_class(self,module,name):
        if (module,name)==('fractions','Fraction'):return Fraction
        raise pickle.UnpicklingError((module,name))


def check():
    name='rational_c4/C0_C2_boundary_exact.pkl'; sha='7c57131cd88155e3fc3e9019332dc5bac9f440d6ebf6af26f586db3907811ca6'
    with zipfile.ZipFile(ROOT/'artifacts/nonperturbative-certificates.zip') as z:data=z.read(name)
    assert hashlib.sha256(data).hexdigest()==sha
    d=RationalWitness(io.BytesIO(data)).load(); T=PureSU2(); w=d['boundary']; c0=d['C0']; c2=d['C2']
    # Independent matrix-basis checks before evaluating the witness.
    assert all(not T.q0(p) for p in T.base.values())
    trace=T.add(T.mul(T.v('ce12',1,0),T.v('ce21',0,1)),T.mul(T.v('ce21',1,0),T.v('ce12',0,1)),T.scale(T.mul(T.v('ch1',1,0),T.v('ch1',0,1)),2)); assert trace==c0
    assert not T.q0(c0) and not T.q0(c2[0]); terms=[]
    for r in range(3):
        for word in w[r]:
            b=sum(u[0].startswith('b') for u in word); assert (b,len(word)-b,b+sum(u[1] for u in word),b+sum(u[2] for u in word))==(2-r,3,4,4)
        assert not T.rotation(w[r],0) and not T.rotation(w[r],1)
        target=T.add(*(star(T,r-i,c0,c2[i]) for i in range(r+1)))
        image=T.add(*(T.qn(r-i,w[i]) for i in range(r+1)))
        assert target==d['products'][r], ('product',r,len(T.add(target,T.scale(d['products'][r],-1))))
        assert target==image, ('differential',r,len(T.add(target,T.scale(image,-1))))
        terms.append(len(target))
    assert T.q0(w[2])==d['last_target']
    return {'witness_sha256':sha,'product_terms':terms,'primitive_terms':list(map(len,w)),'full_boundary_residual_terms':[0,0,0],'rotation_residuals':0,'bare_ghosts_retained':True,'higher_orders':'zero by b-degree','scope':'C0 star C2 only; does not test hidden square, its nonzero witness, or an instanton matrix element'}


if __name__=='__main__':print(json.dumps(check(),indent=2))
