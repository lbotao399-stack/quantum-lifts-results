"""Check the degree, pole-order and Cartan arithmetic in the R-007 discriminant test."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp

from r007_free_radial import target
from r007_soft_form_factor import four_scalar

ROOT=Path(__file__).resolve().parent


def partitions(n, start=1):
    if not n:
        yield ()
        return
    for k in range(start,n+1):
        for tail in partitions(n-k,k):
            yield (k,)+tail


def run():
    source_bytes=(ROOT/'dmk_o330_hw.json').read_bytes()
    assert hashlib.sha256(source_bytes).hexdigest()=='8bee15e4f741312ca2c05e711b6aec8afcd69fe47dddf78e7671a05f5bab9f34'
    source=json.loads(source_bytes)['terms']
    occupations=Counter()
    for coeff,multi in source:
        fields=[l for trace in multi for l in trace]
        degrees=[sum(l==i if l<3 else l-3!=i for l in fields) for i in range(3)]
        assert Fraction(coeff) and len(fields)==10 and degrees==[7,4,4]
        occupations[sum(l!=0 for l in fields)]+=1
    assert len(source)==6139
    # Explicit Lagrange interpolation over the generic traceless SU(3) Cartan.
    x,y,t=sp.symbols('x y t')
    z=(x,y,-x-y)
    vand=sp.Matrix([[a**i for a in z] for i in range(3)])
    inv=vand.inv()
    for i,a in enumerate(z):
        denom=sp.prod(a-b for j,b in enumerate(z) if i!=j)
        coeff=sp.Poly(sp.prod(t-b for j,b in enumerate(z) if i!=j),t)
        assert all(sp.cancel(inv[i,k]*denom-coeff.nth(k))==0 for k in range(3))
        assert all(sp.cancel(coeff.eval(b)/denom-int(i==j))==0 for j,b in enumerate(z))
    delta=sp.prod((a-b)**2 for i,a in enumerate(z) for b in z[i+1:])
    assert sp.expand(delta-(sum(a*a for a in z)**3/2-3*sum(a**3 for a in z)**2))==0
    # KVM (105)-(107): each proper splitting has at least two positive-weight letters.
    poles=[0]
    for w in range(1,9):
        proper=[p for p in partitions(w) if len(p)>=2]
        poles.append(1+max([0]+[sum(poles[i] for i in p) for p in proper]))
        assert poles[w]==2*w-1
    w=8
    max_remainder=max(sum(poles[i] for i in p) for p in partitions(w))
    # At fixed positive-letter count t, h(δh)^k uses at most w-t+1 inversions.
    # A torus-neutral monomial of t letters needs at most t spectral projectors.
    primitive_bound=max((w-s+1)+t for s in range(1,w+1) for t in range(s,w+1))
    assert (max_remainder,primitive_bound)==(15,16)
    clearing=(max_remainder+1)//2+(primitive_bound+1)//2
    assert clearing==16
    steps=[]
    for m in range(5):
        degree=(7-m,4-m,4-m,m,m)
        assert sum(degree[i] for i in (1,2,3,4))==8
        steps.append(dict(order=m,undressed_weight=degree,power=6 if m==0 else clearing))
    exponent=sum(row['power'] for row in steps)
    assert exponent==70
    # All 5,080 terms enter before the Cartan specialization. The four
    # distinguishable external legs give the factor 2!2! in this coefficient.
    words=target()
    assert len(words)==5080
    z0,B,C=(1,2,-3),(1,0,-1),(0,1,-1)
    coefficients,_=four_scalar(words,z0,B,C)
    val=coefficients[(1,1,1,1)]/4
    discr=int(delta.subs({x:z0[0],y:z0[1]}))
    assert val==Fraction(-14112,5) and discr==400
    return dict(source_terms=len(source),positive_letters={str(k):v for k,v in sorted(occupations.items())},
                max_poles_single_diagonal_letter=poles[1:],weight=w,
                graviton_numerator_gap_poles=max_remainder,primitive_gap_poles_after_numerator=primitive_bound,
                general_clearing_exponent=clearing,steps=steps,total_exponent=exponent,
                target_terms=len(words),probe_discriminant=discr,probe_value=str(val),
                dressed_target_value=str(discr**exponent*val))

if __name__=='__main__':
    if not __debug__:
        raise ValueError('enable assertions')
    print(json.dumps(run(),indent=2))
