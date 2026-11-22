"""Exact free-norm projection on a gauge-invariant classical boundary."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction
from functools import lru_cache
import json
from math import prod
from pathlib import Path

import r007_classical_adjoint as adj
import r007_free_radial as radial
import r007_full_radial as full
import r007_primary_replay as jets


def permutation(cycles,length):
    a=[-1]*length
    for tr in cycles:
        for i,j in zip(tr,tr[1:]+tr[:1]):
            a[i]=j
    assert sorted(a)==list(range(length))
    return tuple(a)


def groups(state,keep=None):
    out=defaultdict(Counter)
    for w,c in state.items():
        flat=tuple(x for tr in w for x in tr)
        profile=tuple(sorted(flat))
        if keep is not None and profile not in keep:
            continue
        order=sorted(range(len(flat)),key=flat.__getitem__)
        slot={p:i for i,p in enumerate(order)}
        off=0
        cycles=[]
        for tr in w:
            cycles.append([slot[i] for i in range(off,off+len(tr))])
            off+=len(tr)
        out[profile][permutation(cycles,len(flat))]+=c*full.permutation_sign(order,tuple(map(jets.odd,flat)))
    return out


@lru_cache(100000)
def components(a,b):
    # Simultaneously relabel the two oriented trace permutations. Each
    # connected component is a separate colour sum, including singlet removal.
    untouched=set(range(len(a)))
    ans=[]
    while untouched:
        root=min(untouched)
        comp={root}
        todo=[root]
        for i in todo:
            for v in (a[i],b[i]):
                if v not in comp:
                    comp.add(v)
                    todo.append(v)
        untouched-=comp
        candidates=[]
        for first in comp:
            order=[first]
            ren={first:0}
            for i in order:
                for v in (a[i],b[i]):
                    if v not in ren:
                        ren[v]=len(order)
                        order.append(v)
            candidates.append((tuple(ren[a[i]] for i in order),tuple(ren[b[i]] for i in order)))
        ans.append(min(candidates))
    return tuple(sorted(ans))


@lru_cache(None)
def cycle_masks(a):
    seen=set()
    masks=[]
    for start in range(len(a)):
        if start in seen:
            continue
        x=start
        m=0
        while x not in seen:
            seen.add(x)
            m|=1<<x
            x=a[x]
        masks.append(m)
    return masks


@lru_cache(None)
def reduced(a,mask):
    vals=[]
    for i in range(len(a)):
        if mask>>i&1:
            j=a[i]
            while not mask>>j&1:
                j=a[j]
            vals.append(j)
    return vals,sum(not (x&mask) for x in cycle_masks(a))


@lru_cache(None)
def colour(a,b,n):
    # For a deleted propagator the singlet contributes -1/n. A trace with
    # all its slots deleted contributes n. The U(n) contraction on the
    # retained slots counts cycles of a b^{-1}.
    l=len(a)
    total=0
    for mask in range(1<<l):
        aa,empty1=reduced(a,mask)
        bb,empty2=reduced(b,mask)
        perm=dict(zip(bb,aa))
        loops=0
        for v in aa:
            if v in perm:
                loops+=1
                while v in perm:
                    v=perm.pop(v)
        k=mask.bit_count()
        total+=(-1)**(l-k)*n**(k+empty1+empty2+loops)
    return Fraction(total,n**l)


def inner(a,b,n):
    left=groups(a)
    right=groups(b,keep=set(left))
    total=Fraction()
    for profile in sorted(left.keys()&right.keys()):
        topologies=Counter()
        for m,s in full.matchings(profile):
            for bet,d in right[profile].items():
                relabelled=[-1]*len(m)
                for i,j in enumerate(bet):
                    relabelled[m[i]]=m[j]
                bet_new=tuple(relabelled)
                for al,c in left[profile].items():
                    topologies[components(al,bet_new)]+=s*c*d
        weight=prod(map(radial.inner_letter,profile))
        total+=weight*sum(c*prod(colour(al,be,n) for al,be in comp)
                          for comp,c in topologies.items() if c)
    return total


def run():
    target=radial.target()
    profile=tuple([(0,0,(0,))]*5+[(0,0,(1,)),(0,2,(1,2)),(2,0,(2,))])
    ka={w:c for w,c in adj.adjoint(target,3).items()
        if tuple(sorted(l for tr in w for l in tr))==profile}
    assert len(ka)==6
    qb=jets.act(ka,('Q',0,0))
    assert len(qb)==82 and not jets.act(qb,('Q',0,0)) and not jets.rho(qb)
    knorm=inner(ka,ka,3)
    boundary=inner(qb,qb,3)
    mixed=inner(qb,target,3)
    assert knorm==mixed and min(knorm,boundary)>0
    upper=Fraction(json.loads(Path(__file__).with_name('r007-full-radial.json').read_bytes())['full_SU3_free_norm'])
    lower=Fraction(json.loads(Path(__file__).with_name('r007-free-radial.json').read_bytes())['Cartan_Gram_defining_trace_unit_singleton'])
    projected=upper-mixed*mixed/boundary
    assert lower<=projected<upper
    return dict(primitive_terms=len(ka),boundary_terms=len(qb),primitive_trace_recipe=[[w,str(c)] for w,c in sorted(ka.items())],
                primitive_norm=str(knorm),boundary_norm=str(boundary),boundary_overlap=str(mixed),
                boundary_subtraction_coefficient=str(mixed/boundary),squared_norm_after_one_projection=str(projected),
                lower_bound_harmonic_squared_norm=str(lower),colour_contraction='defining trace, traceless projector at every matched slot')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        raise RuntimeError('run with assertions enabled')
    data=run()
    a.out.write_text(json.dumps(data,indent=2)+'\n')
    print({k:v for k,v in data.items() if k!='primitive_trace_recipe'},flush=True)
