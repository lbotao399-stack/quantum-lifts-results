"""Exact Q3 input-jet moments from pinned currents3 and lowmeasures.

The integral weight is t1^2 t2 t3 t4(1-t4)t6^2 on four intervals
and the triangle t5+t6 <= 1. Its zero-input value is 1/864.
Relative moments preserve the existing path coefficient 1/108.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import factorial, comb

ONE = {(0,)*6: 1}


def add(*terms):
    result = defaultdict(int)
    for coefficient, polynomial in terms:
        for powers, value in polynomial.items():
            result[powers] += coefficient*value
    return {powers: value for powers, value in result.items() if value}


def multiply(*polynomials):
    result = ONE
    for polynomial in polynomials:
        out = defaultdict(int)
        for first, coefficient in result.items():
            for second, value in polynomial.items():
                out[tuple(x+y for x, y in zip(first, second))] += coefficient*value
        result = {powers: value for powers, value in out.items() if value}
    return result


T = [{tuple(int(i == j) for i in range(6)): 1} for j in range(6)]
NT = [add((1, ONE), (-1, t)) for t in T]
A = (
    add((1, NT[0]), (1, multiply(T[0], NT[1], NT[2])),
        (1, multiply(T[0], NT[1], T[2], NT[3], NT[4]))),
    multiply(T[0], T[1], add((1, NT[2]), (1, multiply(T[2], NT[3], NT[4])))),
    multiply(T[0], T[2], T[3], NT[4]),
    multiply(T[0], T[2], T[4]),
)
B = (multiply(NT[1], NT[3], T[5]), multiply(T[1], NT[3], T[5]),
     multiply(T[3], T[5]), NT[5])


def integral(polynomial):
    result = Fraction(0)
    for (i, j, k, l, m, n), value in polynomial.items():
        result += Fraction(value*factorial(m)*factorial(n+2),
            (i+3)*(j+2)*(k+2)*(l+2)*(l+3)*factorial(m+n+4))
    return result


@lru_cache(None)
def moment(exponents):
    assert len(exponents) == 8 and all(isinstance(n, int) and n >= 0 for n in exponents)
    polynomial = multiply(*(p for p, n in zip(A+B, exponents) for _ in range(n)))
    return 864*integral(polynomial)


@lru_cache(None)
def jet_terms(jets):
    assert len(jets) == 4 and all(len(jet) == 2 and min(jet) >= 0 for jet in jets)
    result = defaultdict(Fraction)
    for split in product(*(product(range(i+1), range(j+1)) for i, j in jets)):
        first = tuple(sum(pair[i] for pair in split) for i in (0, 1))
        second = tuple(sum(jet[i]-pair[i] for jet, pair in zip(jets, split)) for i in (0, 1))
        denominator = 1
        for (i, j), (r, s) in zip(jets, split):
            denominator *= factorial(r)*factorial(s)*factorial(i-r)*factorial(j-s)
        exponents = tuple(sum(pair) for pair in split)+tuple(sum(jet)-sum(pair) for jet, pair in zip(jets, split))
        weight = moment(exponents)/denominator
        for k in range(4):
            orders = (first[0]+3-k, first[1]+k, second[0]+k, second[1]+3-k)
            coefficient = (-1)**k*comb(3, k)*weight
            for order in orders:
                coefficient *= factorial(order)
            result[orders] += coefficient
    return tuple((value, (i, j), (k, l)) for (i, j, k, l), value in result.items() if value)
