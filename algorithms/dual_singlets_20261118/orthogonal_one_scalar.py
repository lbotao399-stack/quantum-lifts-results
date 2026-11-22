"""Supercharge descent for S^(k-1) Tr(lambda_+ phi^1) in orthogonal colours."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod
from random import Random

from b2_schur_complex import centralizer, character, partitions
from mixed_incoming import so_canon
from orthogonal_glueball_even_colours import vanishes_on_planes
from orthogonal_glueball_first_pair import GHOSTS, PHI, PfaffianCoefficient
from orthogonal_glueball_relation_span import determinant, weights
from orthogonal_glueball_sources import missing_adjoint_trace
from orthogonal_glueball_third_pair import encode
from orthogonal_trace_dual import pullback
from terminal_commutators import words


def dual(colours, k):
    sign, w = so_canon((GHOSTS,)*(k-1) + ((GHOSTS[0], PHI[0]),))
    first, a = pullback(colours, {w: sign}, 1)
    second, b = pullback(colours, first, -1)
    return dict(colours=colours, k=k, degrees=[k,k-1,1,0,0],first=a,second=b,
                scores=[[encode(w),str(c)] for w,c in sorted(second.items())])


@lru_cache(None)
def character_weights(k, order, full_weight):
    plus, minus = k-order, k-1-order
    scalar = (order+1, order, order)
    if full_weight:
        scalar_terms = [(sum(v,()), Fraction(1,prod(map(centralizer,v))))
                        for v in product(*(partitions(n) for n in scalar))]
    else:
        scalar_terms = [(z,Fraction(c,centralizer(z))) for z in partitions(sum(scalar))
                        if (c:=character(scalar,z))]
    terms = Counter()
    for x,y in product(partitions(plus),partitions(minus)):
        sign = (-1)**(plus+minus-len(x)-len(y))
        coeff = Fraction(sign, centralizer(x)*centralizer(y))
        for z, s in scalar_terms:
            terms[tuple(sorted(x+y+z))] += coeff*s
    return tuple((key,c) for key,c in terms.items() if c)


def dimension(colours, k, order, full_weight=False):
    value = sum(c*missing_adjoint_trace(z,colours) for z,c in character_weights(k,order,full_weight))
    value = Fraction(value)
    assert value.denominator == 1 and value >= 0
    return int(value)


def relation(k,spec,order=2):
    r = determinant(k,spec)
    if order != 2:
        r.project = PfaffianCoefficient(order,0)
        r.letters += order-2
        r.pf_degree += order-2
    # One extra flavour-1 derivative of the scalar determinant projector.
    # The digit of its old monomial accounts for the derivative factorial.
    bit = 1 << (3*spec['axis'])
    old = r.project.weights
    r.project.weights = {code+bit: c*(1+(code//bit)%8) for code,c in old.items()}
    r.project.divisors = r.project.divisors | {a+bit for a in r.project.divisors}
    return r


@lru_cache(None)
def shapes(colours, k, order=2):
    cost = 2*k+order-1-colours
    options = sorted({min(v,v[::-1]) for j in range(2,cost+2)
                      for seq in product('pms', repeat=j) for v in [''.join(seq)]})
    patterns = []
    def visit(start, remaining, paths, plus, minus, scalar):
        if plus < 0 or minus < 0 or scalar < 0:
            return
        if max(plus,minus) <= 3:
            patterns.append(paths)
        if max(0,plus-3)+max(0,minus-3)>2*remaining:
            return
        for i in range(start,len(options)):
            v = options[i]
            if len(v)-1 <= remaining:
                visit(i,remaining+1-len(v),paths+(v,),plus-v.count('p'),
                      minus-v.count('m'),scalar-v.count('s'))
    visit(0,cost,(),k-order,k-order-1,3*order+1)
    return tuple(sorted(patterns,key=lambda x:(len(x),x)))


def specification(index, colours, k, order=2):
    patterns = shapes(colours,k,order)
    paths = patterns[index%len(patterns)]
    letters = ''.join(paths)
    rng = Random(71345051+1009*colours+65537*index)
    return dict(paths=list(paths), matrices=[weights(rng,4) for _ in paths],
                single=[rng.sample(range(3),count-letters.count(s))
                        for s,count in [('p',k-order),('m',k-order-1)]],
                scalar_directions=[[weights(rng,3) for _ in range(path.count('s'))] for path in paths],
                axis=rng.randrange(3))


def evaluate(spec, colours, k, scores, probes=()):
    r = relation(k,spec)
    assert r.pf_degree > colours
    sites = words((k-2,k-3,3,2,2))
    blind = vanishes_on_planes(spec)
    cartan = [str(r.cartan(w)) if not blind and len(w)>colours//2 else '0' for w in sites]
    value = sum(c*r.coefficient(w) for w,c in scores)
    ordered = [str(r.coefficient(w)) for w in probes]
    r.clear()
    return dict(cartan=cartan,second_score=str(value),ordered=ordered)
