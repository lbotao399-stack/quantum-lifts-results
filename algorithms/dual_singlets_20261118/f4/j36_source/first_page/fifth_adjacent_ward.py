"""Classical terminal mixing at the physical charge of the fifth F4 glueball."""
import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
import json
from math import prod
from pathlib import Path

from flint import fmpq_mat

from colour import Colour, Q, add, apply, compact, items, prepend, wedge
from dotted_filtration import pure_current
from glueball_factor_ward import Anchors, BothLegs, negative_sum, solve_constraints
from intermediate_rank import exchange_trace, tensor_multiplicity
from source_symmetry import centralizer, character, multiplicity, partitions, singlet_trace


def exterior_weights(c, degree):
    groups=defaultdict(list)
    for indices in combinations(range(52),degree):
        mask=wedge(indices)[0]
        groups[c.weight(mask)].append(mask)
    return groups


class ThreeScalars(Anchors):
    def __init__(self,c):
        self.c=c
        self.h,=c.high[4]
        self.low=c.low[4,0]
        self.hweight=c.weight(self.h)
        self.lweight=c.weight(self.low)
        self.triples=exterior_weights(c,3)
        self.pairs=exterior_weights(c,2)
        self.scalar_pairs=self.triples[negative_sum(self.hweight,self.lweight)]
        self.index={m:i for i,m in enumerate(self.scalar_pairs)}
        self.mirror=[(i,Q(v)) for i,v in c.data["longest_weyl"]]

    def reflect(self,mask):
        transformed=[self.mirror[i] for i in items(mask)]
        m,s=wedge([i for i,_ in transformed])
        return m,s*prod(v for _,v in transformed)

    def invariant_rows(self):
        c=self.c
        for r,w in enumerate(c.weights):
            for j,v in enumerate(c.weights):
                if w!=negative_sum(v):
                    continue
                row={}
                for pair,x in c.delta[r].items():
                    t,s=prepend(j,pair)
                    if s: add(row,self.index[t],s*x)
                if compact(row):
                    yield ("scalar",r,j),compact(row)
        lo,s=self.reflect(self.h)
        hi,ss=self.reflect(self.low)
        assert lo==self.low and hi==self.h and s*ss==1
        for t in self.scalar_pairs:
            u,phase=self.reflect(t)
            row={self.index[t]:Q(1)}
            add(row,self.index[u],-s*ss*phase)
            if compact(row):
                yield ("exchange",items(t)),compact(row)
        for root,extreme in [(e,self.low) for e in c.data["raising"]]+[(f,self.h) for f in c.data["lowering"]]:
            vector={extreme:1}
            power=0
            while vector:
                vector=apply(lambda m:c.act(root,m),vector)
                power+=1
            # The same simple-root power kills the scalar test after moving
            # the raising operator past the highest (or lowering past lowest) leg.
            weight=tuple(-power*x for x in c.weights[root])
            for t in self.triples.get(weight,()):
                row={t:1}
                for _ in range(power):
                    row=apply(lambda m:c.act(root,m),row)
                if row:
                    yield ("root_power",root,power,items(t)),{self.index[u]:x for u,x in row.items()}

    def compress(self):
        fixed=solve_constraints(self.invariant_rows(),len(self.index))
        m=fmpq_mat(fixed.pop("rational_rows"))
        rr,rank=m.rref()
        pivots=fixed["pivot_columns"]
        assert rank==len(pivots)
        free=[j for j in range(m.ncols()) if j not in pivots]
        basis=fmpq_mat(m.ncols(),len(free))
        for k,j in enumerate(free):
            basis[j,k]=1
            for i,p in enumerate(pivots):
                basis[p,k]=-rr[i,j]
        assert m*basis==fmpq_mat(m.nrows(),len(free))
        self.basis={t:{j:basis[i,j] for j in range(len(free)) if basis[i,j]} for t,i in self.index.items()}
        fixed["free_scalar_triples"]=[items(self.scalar_pairs[j]) for j in free]
        return fixed,len(free)

    def at_low(self,t):
        return self.basis[t]

    def mixed_rows(self,joint):
        c=self.c
        three_minus=[self.low ^ (1<<j) for j in items(self.low)]
        for v in three_minus:
            for r,w in enumerate(c.weights):
                for pair in self.pairs.get(negative_sum(self.hweight,c.weight(v),w),()):
                    row={}
                    for mark,g in c.delta[r].items():
                        x,y=items(mark)
                        for i,j,phase in ((x,y,1),(y,x,-1)):
                            vv,sv=prepend(i,v)
                            t,st=prepend(j,pair)
                            if sv and st:
                                for k,z in self.evaluate(vv,t).items():
                                    add(row,k,phase*sv*st*g*z)
                    if compact(row):
                        yield ("ghost_scalar",items(v),r,items(pair)),compact(row)
        for i in items(self.h):
            u=self.h ^ (1<<i)
            for v in three_minus:
                for r,w in enumerate(c.weights):
                    for t in self.triples.get(negative_sum(c.weight(u),c.weight(v),w),()):
                        row={}
                        for mark,g in c.delta[r].items():
                            x,y=items(mark)
                            for first,second,phase in ((x,y,1),(y,x,-1)):
                                uu,su=prepend(first,u)
                                vv,sv=prepend(second,v)
                                if su and sv:
                                    for k,z in joint.evaluate(uu,vv,t).items():
                                        add(row,k,phase*su*sv*g*z)
                        if compact(row):
                            yield ("opposite_ghosts",items(u),items(v),r,items(t)),compact(row)


def calculate():
    diagrams=[]
    for p in partitions(11):
        shape=p[::-1]
        n=sum(Fraction(character((4,4),a)*character((3,),b)*character(shape,tuple(sorted(a+b))),
                       centralizer(a)*centralizer(b)) for a in partitions(8) for b in partitions(3))
        if n:
            assert n==1
            diagrams.append((shape,multiplicity(shape,11,singlet_trace)))
    assert {s for s,_ in diagrams}=={(7,4),(6,4,1),(5,4,2),(4,4,3)}
    previous=json.loads(Path(__file__).with_name("glueball_factor_ward.json").read_bytes())
    assert previous["upper_bound"]==0 and previous["taylor_shape"]==[5,4,2]
    other=[]
    b,=pure_current(4,0)
    for first,degree in ((7,0),(6,1)):
        for high in pure_current(first,0):
            k,=pure_current(degree,0)
            n=tensor_multiplicity(high,b,k)
            assert n==0
            other.append(dict(positive=high,negative=b,scalar=k,invariants=n))
    k,=pure_current(3,0)
    count=tensor_multiplicity(b,b,k)
    swap=exchange_trace(b,k)
    solver=ThreeScalars(Colour())
    fixed,dimension=solver.compress()
    assert dimension==(count+swap)//2==17
    joint=BothLegs(solver.c,[solver])
    mixed=solve_constraints(solver.mixed_rows(joint),dimension)
    return dict(physical_charges=dict(E=15,J_L=0,J_R=0,R=[5,5,5]), dotted_and_flavour=[4,4,1,1,1],
                terminal_diagrams=diagrams,raw_terminal_dimension=sum(n for _,n in diagrams),
                zero_harmonic_components=other,ordered_harmonic_multiplicity=count,exchange_trace=swap,
                anchor_constraints=fixed,harmonic_dimension=dimension,mixed_brackets=mixed,
                terminal_classical_upper_bound=dimension-mixed["rank"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    record=calculate()
    a.out.write_text(json.dumps(record,indent=2)+"\n")
    print(record["raw_terminal_dimension"],record["harmonic_dimension"],record["mixed_brackets"]["rank"],
          record["terminal_classical_upper_bound"],flush=True)
