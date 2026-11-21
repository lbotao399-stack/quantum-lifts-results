"""Charge-conjugation-resolved SU(N) adjoint invariants by Schur-Weyl duality."""
from collections import Counter
from functools import lru_cache
from math import comb, factorial

import charge_maps  # Loads the campaign's pinned character runtime.
from exact_dims import character, cycle_data, partitions


@lru_cache(None)
def endomorphism_traces(cycles, size):
    # Permutation of the input letters acts by conjugation on each Specht
    # endomorphism algebra. Reversing every trace transposes that algebra.
    squared = tuple(sorted((m if m % 2 else m // 2
                            for m in cycles for _ in range(1 if m % 2 else 2)), reverse=True))
    ordinary = transposed = 0
    for shape in partitions(sum(cycles)):
        if len(shape) <= size:
            ordinary += character(shape, cycles)**2
            transposed += character(shape, squared)
    return ordinary, transposed


@lru_cache(None)
def adjoint_traces(cycles, size):
    # Remove the scalar identity from each letter. Only invariant subsets
    # of positions contribute to the trace of the given permutation.
    choices = {(): 1}
    for length, count in Counter(cycles).items():
        choices = {old+(length,)*k: coefficient * (-1)**(count-k) * comb(count,k)
                   for old,coefficient in choices.items() for k in range(count+1)}
    ordinary = transposed = 0
    for retained, coefficient in choices.items():
        full, reverse = endomorphism_traces(retained, size)
        ordinary += coefficient * full
        transposed += coefficient * reverse
    # Physical colour conjugation sends an adjoint letter to minus its transpose.
    return ordinary, (-1)**sum(cycles) * transposed


@lru_cache(None)
def parity_parts(bosons, fermions, size, with_identity=False):
    terms = {(): 1}
    denominator = 1
    for counts, odd in ((bosons, False), (fermions, True)):
        for count in counts:
            denominator *= factorial(count)
            following = Counter()
            for previous, coefficient in terms.items():
                for cycle, weight in cycle_data(count, odd):
                    key = tuple(sorted(previous+cycle, reverse=True))
                    following[key] += coefficient*weight
            terms = {cycle: weight for cycle,weight in following.items() if weight}
    ordinary = conjugated = 0
    for cycle, coefficient in terms.items():
        if with_identity:
            full, transposed = endomorphism_traces(cycle, size)
            twisted = (-1)**sum(cycle) * transposed
        else:
            full, twisted = adjoint_traces(cycle, size)
        ordinary += coefficient*full
        conjugated += coefficient*twisted
    both = (ordinary+conjugated, ordinary-conjugated)
    assert all(n >= 0 and n % (2*denominator) == 0 for n in both)
    return tuple(n//(2*denominator) for n in both)
