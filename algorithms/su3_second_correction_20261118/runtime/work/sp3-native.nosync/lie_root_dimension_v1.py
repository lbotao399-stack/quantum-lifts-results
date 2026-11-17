"""Integer Weyl constant terms for adjoint tensor invariants.

No floating torus quadrature or rounding is used. The root coordinates
match the public fortuity groupindex implementation; multiplicities use
the existing symmetric/exterior cycle-index coefficients.
"""
import sp_environment
import sys
sys.path.insert(0, str(sp_environment.WORK/'cluster-next-v2.nosync/work/so8-bmn/src'))
from exact_dims import cycle_data
from collections import defaultdict
from functools import lru_cache
from math import factorial


def roots(group):
    if group == 'G2':
        positive = ((1,0),(0,1),(1,1),(2,1),(1,-1),(1,2))
        return tuple(v for root in positive for v in (root, tuple(-x for x in root))), 2, 12
    if group == 'SU3':
        weights = ((1,0),(0,1),(0,0))
        return tuple(tuple(a[i]-b[i] for i in range(2)) for a in weights for b in weights if a != b), 2, 6
    if group == 'Sp2':
        return ((1,1),(1,-1),(-1,1),(-1,-1),(2,0),(-2,0),(0,2),(0,-2)), 2, 8
    raise NotImplementedError(group)


@lru_cache(None)
def measure(group):
    root_list, rank, order = roots(group)
    result = {(0,)*rank:1}
    for root in root_list:
        out = defaultdict(int, result)
        for weight, coefficient in result.items():
            out[tuple(weight[i]+root[i] for i in range(rank))] -= coefficient
        result = {weight:coefficient for weight,coefficient in out.items() if coefficient}
    assert result[(0,)*rank] == order
    return result


@lru_cache(None)
def invariant_trace(cycles, group):
    root_list, rank, order = roots(group)
    zero = (0,)*rank
    character = {zero:1}
    for length in cycles:
        out = defaultdict(int)
        for weight, coefficient in character.items():
            out[weight] += rank*coefficient
            for root in root_list:
                out[tuple(weight[i]+length*root[i] for i in range(rank))] += coefficient
        character = dict(out)
    value = sum(coefficient*character.get(tuple(-x for x in weight),0)
                for weight,coefficient in measure(group).items())
    assert value % order == 0
    return value//order


@lru_cache(None)
def multiplicity_parts(bosons, fermions, group):
    terms = {():1}
    denominator = 1
    for count, odd in [(n,False) for n in bosons if n]+[(n,True) for n in fermions if n]:
        denominator *= factorial(count)
        out = defaultdict(int)
        for previous, coefficient in terms.items():
            for cycle, weight in cycle_data(count, odd):
                out[tuple(sorted(previous+cycle,reverse=True))] += coefficient*weight
        terms = {cycle:weight for cycle,weight in out.items() if weight}
    value = sum(weight*invariant_trace(cycle,group) for cycle,weight in terms.items())
    assert value % denominator == 0
    value //= denominator
    assert value >= 0
    return value
