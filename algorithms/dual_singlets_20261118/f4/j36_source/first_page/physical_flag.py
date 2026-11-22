"""Check that the scalar flag maps nontrivially to the native dotted and flavour singlet."""
from itertools import permutations
import json

from ghost_trace_quotient import sign


def phase(left,right):
    n = 0
    while left:
        bit = left & -left
        n += (right & (bit-1)).bit_count()
        left ^= bit
    return (-1)**n


def multiply(x,y):
    out = {}
    for u,a in x.items():
        for v,b in y.items():
            if not u&v:
                out[u|v] = out.get(u|v,0) + phase(u,v)*a*b
    return {k:v for k,v in out.items() if v}


def determinant(n,scalar):
    terms = {}
    for cols in permutations(range(n)):
        work = {0:sign(cols)}
        for i,j in enumerate(cols):
            value = {(1<<i)|(1<<(8+j)):1,(1<<(4+i))|(1<<(12+j)):1}
            if i == j and i < 3:
                value[0] = scalar[i]
            work = multiply(work,value)
        for k,v in work.items():
            terms[k] = terms.get(k,0)+v
    return {k:v for k,v in terms.items() if v}


def coefficient(scalar):
    half = multiply(determinant(4,scalar),determinant(3,scalar))
    top = 2**16-1
    return sum(a*half.get(top^k,0)*phase(k,top^k) for k,a in half.items())


def calculate():
    base = coefficient((1,1,1))
    scaled = coefficient((2,3,5))
    assert base and scaled == 900*base
    return dict(column_heights=[4,4,3,3],scalar_highest=[4,4,4,2],
                native_dotted=[4,4],native_flavour=[2,2,2],
                exterior_order=["four alpha","four beta","four plus colours","four minus colours"],
                scalar_diagonal=[[1,1,1],[2,3,5]],coefficients=[base,scaled])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
