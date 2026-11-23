#!/usr/bin/env python3
"""Finite, exact interpolation certificate for universal auxiliary identities.

After division by the common highest-weight polynomial, the determinant
derivative identity is a polynomial of total parameter degree at most 3
over the rational-function field of the nine X variables (Leibniz rule).
Twenty unisolvent parameter nodes therefore prove it identically, rather
than merely sample an unbounded family. Divergences have parameter degree1.
No floating point arithmetic or change to physical Q0 is used.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'scalar'))
import pfaffian_scalar_candidates as p


def derivative(poly, index):
    out = Counter()
    for monomial, coefficient in poly.items():
        if monomial[index]:
            exponents = list(monomial)
            exponents[index] -= 1
            out[tuple(exponents)] += coefficient * monomial[index]
    return {m: c for m, c in out.items() if c}


def determinant_derivative(poly):
    out = Counter()
    for perm in permutations(range(3)):
        value = poly
        for i, j in enumerate(perm):
            value = derivative(value, 3*i+j)
        sign = (-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        for m,c in value.items():
            out[m] += sign*c
    return {m:c for m,c in out.items() if c}


def divergence(poly):
    out = Counter()
    for i in range(3):
        for a,b,value in [(0,0,1),(1,2,Fraction(1,2)),(2,1,Fraction(1,2))]:
            term = derivative(derivative(poly,9+3*i+b),3*i+a)
            for m,c in term.items():
                out[m] += value*c
    return {m:c for m,c in out.items() if c}


def monomial_value(monomial, point):
    a,b,k = point
    i,j,l = monomial
    return a**i*b**j*k**l


def main(output):
    started=time.monotonic()
    nodes=[(a,b,k) for a in range(4) for b in range(4-a) for k in range(4-a-b)]
    assert len(nodes)==20
    determinant=p.determinant([[monomial_value(m,node) for m in nodes] for node in nodes])
    assert determinant
    certificates=[]
    for a,b,k in nodes:
        c=k+1
        base=p.mul(p.power(p.D,c-1),p.mul(p.power(p.AA,a),p.power(p.ZZ,b)))
        target=p.mul(p.D,base)
        factor=c*(c+2*b+1)*(c+2*b+2*a+2)
        assert determinant_derivative(target)==p.scale(base,factor),(a,b,c)
        certificates.append({'a':a,'b':b,'c':c,'factor':factor,'exact_polynomial_terms':len(target)})

    contraction_A=p.add(*(p.mul(p.dot(p.X[0],p.X[i]),p.dot(p.X[0],p.Y[i])) for i in range(3)))
    contraction_Z=p.add(*(p.mul(p.dot(p.Z,p.X[i]),p.dot(p.Z,p.Y[i])) for i in range(3)))
    linear_nodes=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
    divergence_checks=[]
    for which in (1,2,3):
        for i,j,k in linear_nodes:
            a=i+(which==2);b=j+(which==3);c=k+1
            base=p.mul(p.power(p.D,c-1),p.mul(p.power(p.AA,a),p.power(p.ZZ,b)))
            if which==1:
                source=p.mul(base,p.PAIR)
                factor=2*a+4*b+3*c+6
            elif which==2:
                source=p.mul(p.power(p.D,c-1),p.mul(p.power(p.AA,a-1),p.mul(p.power(p.ZZ,b),contraction_A)))
                factor=2*a+2*b+c+4
            else:
                source=p.mul(p.power(p.D,c-1),p.mul(p.power(p.AA,a),p.mul(p.power(p.ZZ,b-1),contraction_Z)))
                factor=c
            assert divergence(source)==p.scale(base,factor),(which,a,b,c)
            divergence_checks.append({'source':which,'a':a,'b':b,'c':c,'factor':factor})
    result={'passed':True,'field':'Q',
       'universal_identity':'det(d/dX) [D^c AA^a ZZ^b] = c(c+2b+1)(c+2b+2a+2) D^(c-1) AA^a ZZ^b',
       'parameter_degree_bound':3,'parameter_nodes':certificates,
       'parameter_evaluation_determinant':str(determinant),
       'why_finite_nodes_prove_identity':'Repeated Leibniz expansion after division by D^(c-1) AA^a ZZ^b has total degree at most3 in a,b,c over Q(X). The20-node evaluation matrix is invertible overQ. All20 values vanish as exact full polynomials in X. Thus all rational-function coefficients vanish identically.',
       'divergence_parameter_degree_bound':1,'divergence_checks':divergence_checks,
       'divergence_unisolvence':'For each of the3sources, the4 shifted affine-simplex nodes have determinant1.',
       'seconds':time.monotonic()-started}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('parameter_nodes','divergence_checks')},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=HERE/'UNIVERSAL_AUXILIARY_IDENTITY.json')
    args=parser.parse_args();main(args.out)
