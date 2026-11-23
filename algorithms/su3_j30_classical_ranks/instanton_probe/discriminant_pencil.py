"""Check the flavour cubic, the discriminant clearing identity and six joint collision probes."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
import json
from math import comb

from cubic_colours import source


@cache
def necklace(word):
    odd = sum(x >= 3 for x in word)
    forms = {}
    sign = 1
    for _ in word:
        if word in forms and forms[word] != sign:
            return (), 0
        forms[word] = sign
        sign *= (-1)**((word[0] >= 3) * (odd - 1))
        word = word[1:] + word[:1]
    least = min(forms)
    return least, forms[least]


def shear(terms, i, j):
    # The elementary flavour matrix acts on scalars and on their dual matter fermions.
    out = Counter()
    for coefficient, monomial in terms:
        value = Fraction(coefficient)
        monomial = tuple(tuple(w) for w in monomial)
        for a, old in enumerate(monomial):
            for b, letter in enumerate(old):
                if letter not in (j, 3+i):
                    continue
                new, sign = (i, 1) if letter == j else (3+j, -1)
                trace, cyclic = necklace(old[:b] + (new,) + old[b+1:])
                if not cyclic:
                    continue
                traces = monomial[:a] + (trace,) + monomial[a+1:]
                odd = [t for t in traces if sum(x >= 3 for x in t) % 2]
                if len(odd) != len(set(odd)):
                    continue
                swaps = sum(u > t for k, u in enumerate(odd) for t in odd[k+1:])
                out[tuple(sorted(traces))] += value * sign * cyclic * (-1)**swaps
    return {m: x for m, x in out.items() if x}


def matchings(labels):
    if not labels:
        yield ()
        return
    for j in range(1, len(labels)):
        for tail in matchings(labels[1:j] + labels[j+1:]):
            yield ((labels[0], labels[j]),) + tail


def det(rows):
    a, b, c = rows
    return (a[0]*(b[1]*c[2] - b[2]*c[1]) + a[1]*(b[2]*c[0] - b[0]*c[2])
            + a[2]*(b[0]*c[1] - b[1]*c[0]))


def verify():
    terms = source()
    weights = {(1,0,0), (0,1,0), (0,0,1)}
    scalar = sorted(weights, reverse=True)
    degree = [*scalar, *(tuple(1-x for x in w) for w in scalar)]
    assert all(tuple(sum(degree[x][i] for t in traces for x in t) for i in range(3)) == (7,4,4)
               for _, traces in terms)
    actions = [(0,1), (1,2), (2,1), (0,2)]
    assert all(not shear(terms,i,j) for i,j in actions)
    coefficients = [(-1)**j * comb(69+j,j) for j in range(4)]
    polynomial = [sum(coefficients[j] * comb(70,k-j) for j in range(k+1)) for k in range(4)]
    assert polynomial == [1,0,0,0]
    directions = (0,1,2,3,4,6)
    assert all(det([(1,t,t*t) for t in triple]) for triple in combinations(directions,3))
    tests = []
    for pairs in matchings(directions):
        lines = [(u*v,-(u+v),1) for u,v in pairs]
        determinant = det(lines)
        assert determinant
        tests.append(dict(matching=pairs, determinant=determinant))
    assert len(tests) == 15
    # Explicit traceless triple of joint eigenvalue vectors. Its nonzero pair differences
    # make all of the first five projected discriminants vanish, while the sixth is nonzero.
    def value(vector,t): return vector[0]+vector[1]*t+vector[2]*t*t
    a,b=(0,-1,1),(6,-5,1)
    # The third secant goes through t=4: multiply a by 1 and b by -6.
    b=tuple(-6*x for x in b)
    eig=[tuple(2*x+y for x,y in zip(a,b)),tuple(-x+y for x,y in zip(a,b)),
         tuple(-x-2*y for x,y in zip(a,b))]
    assert len(set(eig))==3 and all(sum(w[i] for w in eig)==0 for i in range(3))
    discs=[]
    for t in directions:
        x,y,z=[value(v,t) for v in eig];discs.append(((x-y)*(y-z)*(x-z))**2)
    assert discs[:5] == [0]*5 and discs[5]
    return dict(source_terms=len(terms),weight=[7,4,4],shear_residuals=[[i,j,0] for i,j in actions],
                clearing_coefficients=coefficients,truncated_product=polynomial,
                original_sufficient_power=70,uniform_sufficient_power=73,
                directions=directions,matching_tests=tests,five_test_control=dict(eigenvalues=eig,discriminants=discs),
                nonperturbative_lift='open')

if __name__ == '__main__':
    if not __debug__: raise SystemExit('Enable assertions')
    print(json.dumps(verify(),indent=2))
