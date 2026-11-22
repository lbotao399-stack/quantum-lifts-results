"""Finite F4 quartic constraints on the fourteen formal scalar commutator directions."""
from fractions import Fraction
from itertools import permutations, product
import json
from random import Random

import numpy as np
from sympy import Matrix

from flag_ghost_rows import necklaces
from single_trace_flag import adjoint, native26


def multiply(a,b):
    row = {}
    for u,v in a.items():
        for w,y in b.items():
            row[u+w] = row.get(u+w,0)+v*y
    return row


def bracket(a,b):
    out = multiply(a,b)
    for w,y in multiply(b,a).items():
        out[w] = out.get(w,0)-y
    return {w:y for w,y in out.items() if y}


def verify_quartic(rep):
    size = rep.shape[1]
    h = rep[0]
    t2 = int(np.trace(h@h))
    ratio = Fraction(int(np.trace(h@h@h@h)),t2*t2)
    # All coefficients of the degree-four polynomial on the Cartan.
    tests = 0
    for a in range(4):
        for b in range(a,4):
            for c in range(b,4):
                for d in range(c,4):
                    matrices = [rep[j] for j in (a,b,c,d)]
                    left = sum(int(np.trace(matrices[0]@matrices[i]@matrices[j]@matrices[k]))
                               for i,j,k in permutations((1,2,3)))
                    right = sum(int(np.trace(matrices[0]@matrices[i]))*int(np.trace(matrices[j]@matrices[k]))
                                for i,j,k in ((1,2,3),(2,1,3),(3,1,2)))
                    assert left == 2*ratio*right
                    tests += 1
    assert tests == 35
    # A noncommuting scalar triple in the centralizer of a Cartan ghost.
    probe = None
    native = adjoint()
    for i in range(4,52):
        if native[0,i,i]:
            continue
        for j in range(4,52):
            C = rep[i]@rep[j]-rep[j]@rep[i]
            t = int(np.trace(C@C))
            if t and not np.any(h@rep[j]-rep[j]@h):
                value = int(np.trace(h@h@C@C))
                assert value == ratio*t2*t/3
                probe = dict(ghost=0,scalars=[i,j],trace_C_squared=t,ghost_trace=value)
                break
        if probe:
            break
    assert probe
    return dict(representation=size,trace_ratio=str(ratio),polarization=str(ratio/3),
                cartan_coefficients=tests,centralizer_probe=probe)


def relation_row(sample,polarization):
    rational = Fraction(polarization)
    n,d = rational.numerator,rational.denominator
    rng = Random(361400+sample)
    variables = sorted({necklaces(w)[0] for size in range(1,7) for w in product(range(3),repeat=size)
                        if necklaces(w)[1]})
    plain = {key:rng.randint(-3,3) for key in variables}
    dotted = {key:rng.randint(-3,3) for key in variables}
    S = 1+sample%3
    def trace(source,*factors):
        row = {():1}
        for factor in factors:
            row = multiply(row,factor)
        total = 0
        for w,y in row.items():
            tag,s = necklaces(w)
            if s:
                total += y*s*source[tag]
        return total
    X = [{(i,):1} for i in range(3)]
    C = [bracket(X[(i+1)%3],X[(i+2)%3]) for i in range(3)]
    t = trace(plain,C[0],X[0])
    g = trace(dotted,C[0],X[0])
    assert all(trace(plain,C[i],X[j]) == (t if i==j else 0) for i,j in product(range(3),repeat=2))
    a = d*S**2*g*g - n*S**3*t*g
    b = d*S**3*t*g - n*S**4*t*t
    c = e = f = h = 0
    for i,j in product(range(3),repeat=2):
        cc = trace(plain,C[i],C[j])
        xx = trace(plain,X[i],X[j])
        tg = trace(dotted,X[i],X[j])
        c += S**2*(d*trace(dotted,C[i],C[j])-n*S*cc)*tg
        # Odd one-ghost trace generators are specialized to zero only in this
        # independence test. The equality on F4 includes both crossed products.
        e += S**3*cc*(d*tg-n*S*xx)
        terms = [C[i],C[j],X[i],X[j]]
        sym = sum(trace(plain,terms[0],terms[k],terms[l],terms[m])
                  for k,l,m in permutations((1,2,3)))
        f += S**4*(d*sym-6*n*(cc*xx+trace(plain,C[i],X[i])*trace(plain,C[j],X[j])
                                      +trace(plain,C[i],X[j])*trace(plain,C[j],X[i])))
        left,right = bracket(C[i],X[j]),bracket(C[j],X[i])
        h += S**3*(d*trace(dotted,left,right)-n*S*trace(plain,left,right))
    return [a,b,c,e,f,h]


def calculate():
    records = []
    for rep in (native26(),adjoint()):
        r = verify_quartic(rep)
        matrix = Matrix([relation_row(i,r["polarization"]) for i in range(12)])
        basis = list(matrix.T.rref()[1])
        assert len(basis) == 6
        minor = matrix[basis,:]
        r.update(moment_seeds=basis,minor=[[int(x) for x in row] for row in minor.tolist()],
                 determinant=str(minor.det()))
        records.append(r)
    return dict(representations=records,independent_relations_in_each_trace_presentation=6,
                finite_scalar_commutator_rank_upper_in_each=14-6)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
