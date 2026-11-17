"""Exact invariant multiplicities for the adjoint of Sp(rank).

Adjoint = Sym^2(fundamental). The Sp-invariant Schur modules have even
column lengths: partitions with equal adjacent row lengths. The power-sum
plethysm is p_k[h_2] = (p_k^2+p_(2k))/2. All arithmetic is integral.
"""
import sp_environment
import sys
sys.path.insert(0,str(sp_environment.WORK/'cluster-next-v2.nosync/work/so8-bmn/src'))
from exact_dims import partitions,character,cycle_data
from collections import defaultdict
from functools import lru_cache
from math import factorial


@lru_cache(None)
def invariant_trace(cycles,rank):
    degree = sum(cycles)
    shapes = [tuple(value for part in shape for value in (part,part))
              for shape in partitions(degree) if len(shape)<=rank]
    terms = {():1}
    for length in cycles:
        new = defaultdict(int)
        for old,coefficient in terms.items():
            new[tuple(sorted(old+(length,length),reverse=True))] += coefficient
            new[tuple(sorted(old+(2*length,),reverse=True))] += coefficient
        terms = dict(new)
    value = sum(coefficient*sum(character(shape,cycles0) for shape in shapes)
                for cycles0,coefficient in terms.items())
    denominator = 2**len(cycles)
    assert value % denominator == 0
    return value//denominator


@lru_cache(None)
def multiplicity_parts(bosons,fermions,rank):
    terms = {():1}
    denominator = 1
    for count,odd in [(count,False) for count in bosons]+[(count,True) for count in fermions]:
        denominator *= factorial(count)
        new = defaultdict(int)
        for old,coefficient in terms.items():
            for cycle,weight in cycle_data(count,odd):
                new[tuple(sorted(old+cycle,reverse=True))] += coefficient*weight
        terms = {cycle:weight for cycle,weight in new.items() if weight}
    value = sum(weight*invariant_trace(cycle,rank) for cycle,weight in terms.items())
    assert value % denominator == 0
    value //= denominator
    assert value>=0
    return value


def multiplicity(content,rank):
    bosons = tuple(sorted(content[:3]+(content[6],),reverse=True))
    fermions = tuple(sorted(content[3:6],reverse=True))
    return multiplicity_parts(bosons,fermions,rank)
