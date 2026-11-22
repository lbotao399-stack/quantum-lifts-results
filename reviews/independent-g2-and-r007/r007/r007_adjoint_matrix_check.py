"""Reproduce one boundary overlap by explicit traceless 3-by-3 generators."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from random import Random

import r007_free_radial as radial
import r007_primary_replay as jets


def multiply(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))


@lru_cache(40000)
def matrix_trace(mats):
    m=mats[0]
    for t in mats[1:]:
        m=multiply(m,t)
    return sum(m[i][i] for i in range(3))


def entry(i,j):
    return tuple(tuple(int(k==i and l==j) for l in range(3)) for k in range(3))


def diagonal(x):
    return tuple(tuple(x[i] if i==j else 0 for j in range(3)) for i in range(3))


def trace_basis():
    basis=[(entry(i,j),entry(j,i),1) for i in range(3) for j in range(3) if i!=j]
    for vals,n in (((1,-1,0),2),((1,1,-2),6)):
        h=diagonal(vals)
        basis.append((h,h,n))
    assert all(matrix_trace((t,s))/F(n)==int(i==j)
               for i,(t,_,_) in enumerate(basis) for j,(_,s,n) in enumerate(basis))
    return basis


def rand_matrices(profile,seed):
    # The last species is odd, so its polarization uses one Grassmann marker.
    rng=Random(seed)
    out={}
    for parity in (0,1):
        for key in sorted(x for x in set(profile) if jets.odd(x)==parity):
            m=[[rng.randrange(-2,3) for _ in range(3)] for _ in range(3)]
            m[-1][-1]=-m[0][0]-m[1][1]
            out[key]=tuple(tuple(row) for row in m)
    return out


def profile_value(state,profile,seed):
    # Pair the nine-letter target with Q0 applied to an eight-letter coherent
    # test coefficient of these matrices. Every removed pair is summed in
    # the eight-generator trace-dual basis, without a trace sewing identity.
    assert sum(map(jets.odd,profile))==1
    matrices=rand_matrices(profile,seed)
    census=Counter(profile)
    basis=trace_basis()
    total=F()
    routed=0
    for w,c in state.items():
        letters=tuple(x for tr in w for x in tr)
        lengths=tuple(map(len,w))
        for i,u in enumerate(letters):
            for j,v in enumerate(letters):
                if i==j or set(u[2])&set(v[2]):
                    continue
                joined=(u[0]+v[0],u[1]+v[1],tuple(sorted(u[2]+v[2])))
                if joined not in census:
                    continue
                remaining=[x for k,x in enumerate(letters) if k!=i and k!=j]
                if Counter([joined]+remaining)!=census:
                    continue
                split=sum(k for k,split in jets.q_letter(joined) if split==(u,v))
                if not split:
                    continue
                routed+=1
                h=radial.inner_letter(u)*radial.inner_letter(v)/radial.inner_letter(joined)
                # Moving the ordered pair to the front gives the Wick sign.
                order=[i,j]+[k for k in range(9) if k not in (i,j)]
                odd=[k for k in order if jets.odd(letters[k])]
                sign=(-1)**sum(a>b for k,a in enumerate(odd) for b in odd[k+1:])
                colour=F()
                for tu,dual_u,nu in basis:
                    for tv,dual_v,nv in basis:
                        vtx=matrix_trace((matrices[joined],dual_v,dual_u))
                        if not vtx:
                            continue
                        slots=tuple(tu if k==i else tv if k==j else matrices[x] for k,x in enumerate(letters))
                        value=vtx
                        start=0
                        for size in lengths:
                            value*=matrix_trace(slots[start:start+size])
                            start+=size
                        colour+=F(value,nu*nv)
                total+=c*split*h*sign*colour
    return total,routed


def run():
    profile=(((0,0,(0,)),)*5+
             ((0,0,(1,)),(0,2,(1,2)),(2,0,(2,))))
    val,routed=profile_value(radial.target(),profile,596023)
    assert val==912
    return dict(seed=596023,traceless_colour_basis=8,matched_ordered_splits=routed,
                profile=profile,direct_boundary_overlap=str(val))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        raise RuntimeError('run with assertions enabled')
    data=run()
    a.out.write_text(json.dumps(data,indent=2)+'\n')
    print(data,flush=True)
