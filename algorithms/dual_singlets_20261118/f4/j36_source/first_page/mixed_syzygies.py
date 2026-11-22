"""Resolve cancellations between the three matter-fermion flavours at low scalar degrees."""
import argparse
from array import array
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations_with_replacement, product
from math import lcm
import json
import os
import subprocess
import sys

from colour import Colour, Q, add, compact, items
from one_fermion import orbits

PAIRS = [tuple(k for k in range(3) if k!=i) for i in range(3)]
PRIME = 1000003

def remove(degrees,labels):
    left = list(degrees)
    for i in labels:
        left[i] -= 1
    return tuple(left) if min(left)>=0 else None

def with_scalars(spectators,*extra):
    return tuple(sorted((*spectators,*extra)))

class Relations:
    def __init__(self):
        self.c = Colour()
        denominator = lcm(*(int(Q(x).denominator) for row in self.c.delta for x in row.values()),
                          *(int(x.denominator) for x in self.c.inverse.values()))
        self.delta = [{m:int(x*denominator) for m,x in row.items()} for row in self.c.delta]
        self.inverse = {k:int(x*denominator) for k,x in self.c.inverse.items()}
        self.by_weight = defaultdict(list)
        for i,w in enumerate(self.c.weights):
            self.by_weight[w].append(i)

    def assignments(self,slots,weight):
        """Colour labels at fixed weight, sorted only within each scalar flavour."""
        if not slots:
            if not any(weight):
                yield ()
            return
        for head in product(range(52),repeat=len(slots)-1):
            if any(slots[i]==slots[i+1] and head[i]>head[i+1] for i in range(len(head)-1)):
                continue
            needed = tuple(k-sum(self.c.weights[a][j] for a in head) for j,k in enumerate(weight))
            for a in self.by_weight.get(needed,()):
                if head and slots[-1]==slots[-2] and head[-1]>a:
                    continue
                yield (*head,a)

    @cache
    def polarized(self,colour,flavours):
        terms=[]
        for mask,z in self.delta[colour].items():
            a,b = items(mask)
            i,j = flavours
            terms.append((((i,a),(j,b)),z))
            terms.append((((i,b),(j,a)),-z))
        return terms

    @staticmethod
    def key(spect):
        # Each flavour has at most two scalar slots. Their sum and sum of
        # squares determine the unordered colours, including colour zero.
        return sum((a+128*a*a) << (20*i) for i,a in spect)

    @cache
    def tree(self,a,flavours):
        return [(self.key(extra),z) for extra,z in self.polarized(a,flavours)]

    def block(self,degrees,w):
        def fields(counts):
            return tuple(i for i,n in enumerate(counts) for _ in range(n))
        cols = []
        for i,pair in enumerate(PAIRS):
            remaining = remove(degrees,pair)
            if remaining is None:
                continue
            slots = fields(remaining)
            for a,*spect in self.assignments((-1,*slots),w):
                cols.append((i,a,tuple(zip(slots,spect))))
        idx={c:j for j,c in enumerate(cols)}
        images = []
        for i,a,spect in cols:
            row = {}
            fixed = self.key(spect)
            for extra,z in self.tree(a,PAIRS[i]):
                add(row,fixed+extra,z)
            images.append(compact(row))
        boundaries,products = [],[]
        # Two-fermion and three-theta jets furnish the classical predecessors.
        for i,j in combinations_with_replacement(range(3),2):
            rem = remove(degrees,PAIRS[i]+PAIRS[j])
            if rem is None:
                continue
            slots = fields(rem)
            for a,b,*v in self.assignments((-2,-1,*slots),w):
                if i==j and a>=b:
                    continue
                spect = tuple(zip(slots,v))
                row = {}
                for k,x,l,y,sg in ((i,a,j,b,1),(j,b,i,a,-1)):
                    for extra,z in self.polarized(x,PAIRS[k]):
                        add(row,idx[l,y,with_scalars(spect,*extra)],sg*z)
                boundaries.append(compact(row))
        rem = remove(degrees,(0,1,2))
        if rem is not None:
            slots = fields(rem)
            for a,*v in self.assignments((-1,*slots),w):
                spect = tuple(zip(slots,v))
                row = {}
                for mask,z in self.delta[a].items():
                    b,k = items(mask)
                    for i,sg in enumerate((1,-1,1)):
                        for r,s,sign in ((b,k,sg),(k,b,-sg)):
                            add(row,idx[i,r,with_scalars(spect,(i,s))],z*sign)
                boundaries.append(compact(row))
        # Every component below belongs to the classical traceless bilinear.
        grav = [[(i,j,1)] for i,pair in enumerate(PAIRS) for j in pair]
        grav += [[(0,0,1),(1,1,1)],[(1,1,1),(2,2,1)]]
        for pattern in grav:
            i,j,sg = pattern[0]
            rem = remove(degrees,PAIRS[i]+(j,))
            if rem is None:
                continue
            slots = fields(rem)
            for v in self.assignments(slots,w):
                spect = tuple(zip(slots,v))
                row={}
                for i,j,sg in pattern:
                    for (a,b),z in self.inverse.items():
                        add(row,idx[i,a,with_scalars(spect,(j,b))],z*sg)
                products.append(compact(row))
        for vec in boundaries+products:
            check={}
            for j,x in vec.items():
                for s,y in images[j].items():
                    add(check,s,x*y)
            assert not compact(check)
        tree_rank=rank_mod(images)
        consecutive=rank_mod(boundaries)
        joint=rank_mod(boundaries+products)
        assert tree_rank+joint<=len(cols)
        return dict(columns=len(cols),classical_rank_lower=tree_rank,boundary_rank_lower=consecutive,
                    known_cycle_rank_lower=joint,unresolved_upper=len(cols)-tree_rank-joint)

@cache
def residue(x):
    if isinstance(x,int):
        return x%PRIME
    x=Q(x)
    return int(x.numerator)*pow(int(x.denominator),-1,PRIME)%PRIME

def echelon_rank(rows):
    freq=Counter(t for row in rows for t in row)
    order={k:i for i,k in enumerate(sorted(freq,key=lambda k:(freq[k],k)))}
    pivots={}
    for input_row in rows:
        row={order[t]:residue(x) for t,x in input_row.items() if residue(x)}
        while row:
            k=min(row)
            x=row[k]
            if k not in pivots:
                inv=pow(x,-1,PRIME)
                pivots[k]={t:(v*inv)%PRIME for t,v in row.items()}
                break
            for t,v in pivots[k].items():
                y=(row.get(t,0)-x*v)%PRIME
                if y:
                    row[t]=y
                else:
                    row.pop(t,None)
    return len(pivots)

def rank_mod(rows):
    """Eliminate output coordinates incident on one or two input columns before echelon reduction."""
    program = os.environ.get("N4_SPARSE_RANK")
    if program:
        coordinates = {}
        data = array("I",[len(rows),0])
        for row in rows:
            offset = len(data)
            data.append(0)
            for key,x in row.items():
                v=residue(x)
                if v:
                    data.extend((coordinates.setdefault(key,len(coordinates)),v))
            data[offset]=(len(data)-offset-1)//2
        data[1]=len(coordinates)
        assert data.itemsize==4
        if sys.byteorder!="little":
            data.byteswap()
        return int(subprocess.check_output([program],input=data.tobytes()))
    rank = 0
    while rows:
        frequency = Counter(k for row in rows for k in row)
        short = {k:[] for k,n in frequency.items() if n<=2}
        if not short:
            break
        for i,row in enumerate(rows):
            for k,x in row.items():
                if k in short:
                    v = residue(x)
                    if v:
                        short[k].append((i,v))
        # Each variable is its coefficient times its root variable. A forced
        # zero propagates to the whole connected component of two-term equations.
        parent = list(range(len(rows)))
        coefficient = [1]*len(rows)
        size = [1]*len(rows)
        zero = set()
        def find(i):
            if parent[i]!=i:
                root,c = find(parent[i])
                coefficient[i] = coefficient[i]*c%PRIME
                parent[i] = root
            return parent[i],coefficient[i]
        for equation in short.values():
            if not equation:
                continue
            (i,a),*rest = equation
            u,c = find(i)
            if not rest:
                zero.add(u)
                continue
            j,b = rest[0]
            v,d = find(j)
            a,b = a*c%PRIME,b*d%PRIME
            if u==v:
                if (a+b)%PRIME:
                    zero.add(u)
                continue
            if size[u]>size[v]:
                u,v,a,b = v,u,b,a
            parent[u] = v
            coefficient[u] = -b*pow(a,-1,PRIME)%PRIME
            size[v] += size[u]
            if u in zero:
                zero.add(v)
        alive = {i for i in range(len(rows)) if parent[i]==i and i not in zero}
        if len(alive)==len(rows):
            break
        rank += len(rows)-len(alive)
        combined = {i:{} for i in sorted(alive)}
        for i,row in enumerate(rows):
            u,c = find(i)
            if u not in alive:
                continue
            out = combined[u]
            for k,x in row.items():
                if k not in short:
                    out[k] = (out.get(k,0)+c*residue(x))%PRIME
        rows = [compact(row) for row in combined.values() if row]
    return rank + echelon_rank(rows)

def verify(degrees, only=None):
    engine=Relations()
    rows=[]
    for w,size in orbits(engine.c,sum(degrees)-1):
        if only is not None and w!=only:
            continue
        block=engine.block(degrees,w)
        if block['columns']:
            rows.append(dict(weight=w,orbit=size,**block))
    keys=('columns','classical_rank_lower','boundary_rank_lower','known_cycle_rank_lower','unresolved_upper')
    out={k:sum(b['orbit']*b[k] for b in rows) for k in keys}
    expected={(1,1,1):3*52**2,(2,1,1):52*(52*53//2)+2*52**3,
              (2,2,1):52**4+2*52**2*(52*53//2)}
    if only is None:
        assert out['columns']==expected[degrees]
    else:
        assert len(rows)==1
    return dict(flavour_degrees=degrees,prime=PRIME,blocks=sorted(rows,key=lambda b:b['weight']),total=out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('degrees',choices=('111','211','221'))
    p.add_argument('--weight',type=int,nargs=4)
    args=p.parse_args()
    if not __debug__:
        p.error('Run with assertions enabled')
    selected = tuple(args.weight) if args.weight is not None else None
    print(json.dumps(verify(tuple(map(int,args.degrees)),selected),indent=2))
