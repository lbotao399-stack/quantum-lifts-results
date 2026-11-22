"""Resolve cancellations between the three matter-fermion flavours at low scalar degrees."""
import argparse
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, combinations_with_replacement, product
import json

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

    def polarized(self,colour,flavours):
        for mask,z in self.c.delta[colour].items():
            a,b = items(mask)
            i,j = flavours
            yield ((i,a),(j,b)),z
            yield ((i,b),(j,a)),-z

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
            for extra,z in self.polarized(a,PAIRS[i]):
                add(row,with_scalars(spect,*extra),z)
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
                for mask,z in self.c.delta[a].items():
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
                    for (a,b),z in self.c.inverse.items():
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
    x=Q(x)
    return int(x.numerator)*pow(int(x.denominator),-1,PRIME)%PRIME

def rank_mod(rows):
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

def verify(degrees):
    engine=Relations()
    rows=[]
    # Two colour slots for total degree three, three for total degree four.
    # All adjoint tensor-square weights also occur in the tensor cube (zero weight).
    for w,size in orbits(engine.c):
        block=engine.block(degrees,w)
        if block['columns']:
            rows.append(dict(weight=w,orbit=size,**block))
    keys=('columns','classical_rank_lower','boundary_rank_lower','known_cycle_rank_lower','unresolved_upper')
    out={k:sum(b['orbit']*b[k] for b in rows) for k in keys}
    expected={(1,1,1):3*52**2,(2,1,1):52*(52*53//2)+2*52**3}
    assert out['columns']==expected[degrees]
    return dict(flavour_degrees=degrees,prime=PRIME,blocks=sorted(rows,key=lambda b:b['weight']),total=out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('degrees',choices=('111','211'))
    args=p.parse_args()
    if not __debug__:
        p.error('Run with assertions enabled')
    print(json.dumps(verify(tuple(map(int,args.degrees))),indent=2))
