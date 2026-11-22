"""Convert polarized trace moments through degree eight to the native 26-dimensional traces."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
import json
from math import factorial, prod

import numpy as np
from sympy import Matrix

from single_trace_flag import adjoint, native26


def polynomial(weights, degree):
    out = {}
    for places in combinations_with_replacement(range(4),degree):
        powers = tuple(places.count(i) for i in range(4))
        value = sum(prod(v**p for v,p in zip(weight,powers)) for weight in weights)
        if value:
            out[powers] = value*factorial(degree)//prod(factorial(p) for p in powers)
    return out


def multiply(left,right):
    out = Counter()
    for p,x in left.items():
        for q,y in right.items():
            out[tuple(a+b for a,b in zip(p,q))] += x*y
    return {k:x for k,x in out.items() if x}


def cartan(rep):
    vectors = [r.diagonal() for r in rep[:4]]
    assert all(np.array_equal(np.diag(w),r) for w,r in zip(vectors,rep[:4]))
    return tuple(tuple(map(int,v)) for v in zip(*vectors))


def calculate():
    minimal, adj = native26(), adjoint()
    weights = cartan(minimal)
    adj_weights = cartan(adj)
    wedge_weights = Counter(tuple(x+y for x,y in zip(a,b)) for a,b in combinations(weights,2))
    wedge_weights.subtract(adj_weights)
    assert min(wedge_weights.values()) >= 0 and sum(wedge_weights.values()) == 273
    tests = []
    for rep_weights in (adj_weights,tuple(wedge_weights.elements())):
        for n in (1,3,5,7):
            assert not polynomial(rep_weights,n) and not polynomial(weights,n)
        outputs = {}
        for n in (2,4,6,8):
            parts = [w for k in range(1,5) for w in combinations_with_replacement((2,6,8),k) if sum(w)==n]
            columns = []
            for pattern in parts:
                value = {(0,0,0,0):1}
                for deg in pattern:
                    value = multiply(value,polynomial(weights,deg))
                columns.append(value)
            row = polynomial(rep_weights,n)
            keys = sorted(set(row).union(*(col.keys() for col in columns)))
            a = Matrix([[column.get(k,0) for column in columns] for k in keys])
            b = Matrix([row.get(k,0) for k in keys])
            chosen = list(a.T.rref()[1])
            minor = a[chosen,:]
            assert len(chosen) == len(parts) and minor.det()
            coefficients = minor.inv()*b[chosen,:]
            assert a*coefficients == b
            outputs[n] = dict(monomials=parts, coefficients=list(map(str,coefficients)),
                              independent_cartan_monomials=[keys[i] for i in chosen],
                              minor=str(minor.det()))
        tests.append(dict(dimension=len(rep_weights),degrees=outputs))

    # Polarization on pairs of noncommuting native generators, coefficient by
    # coefficient in two variables, independently of the Cartan expansion.
    for i,j in ((4,51),(13,46),(19,50),(21,44)):
        x = np.array([minimal[i]+minimal[0],minimal[j]+minimal[1]],dtype=object)
        y = np.array([adj[i]+adj[0],adj[j]+adj[1]],dtype=object)
        def traces(pair):
            current = [np.eye(pair.shape[-1],dtype=object)]
            out = {}
            for n in range(1,9):
                nxt = [m@pair[0] for m in current]+[np.zeros_like(current[0])]
                for k,m in enumerate(current):
                    nxt[k+1] += m@pair[1]
                current = nxt
                out[n] = [int(m.trace()) for m in current]
            return out
        left,right = traces(x),traces(y)
        for n,data in tests[0]["degrees"].items():
            combined = [0]*(n+1)
            for pattern,c in zip(data["monomials"],data["coefficients"]):
                term = [1]
                for degree in pattern:
                    term = np.convolve(np.array(term,dtype=object),left[degree]).tolist()
                combined = [a+Fraction(c)*b for a,b in zip(combined,term)]
            assert combined == right[n]
    return dict(primitive_degrees_through_eight=[2,6,8],representations=tests,
                native_two_variable_checks=4,tested_degrees=[2,4,6,8])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
