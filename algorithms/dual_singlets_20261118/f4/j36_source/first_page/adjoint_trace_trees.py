"""Eliminate an adjoint-contracted leaf of a defining-trace tree in F4."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import permutations
import json
from math import factorial, lcm

from flint import fmpz_mat
import numpy as np

from single_trace_flag import native26
from source_symmetry import character


def trace_pair(a, b):
    product = a * b
    return sum(product[i, i] for i in range(product.nrows()))


@cache
def projection(slots):
    """Return (Lie word, trace pair) terms in the trace-orthogonal projection."""
    n = len(slots)
    out = Counter()
    if n == 1:
        return {(slots[0], ()): Fraction(1)}
    if n == 3:
        for i in range(3):
            out[slots[i], slots[:i] + slots[i+1:]] += Fraction(1, 36)
    # Subtract each permutation from the original order by adjacent exchanges.
    # Their average leaves the fully symmetric part, handled above.
    for permutation in permutations(range(n)):
        order = list(range(n))
        for i, wanted in enumerate(permutation):
            j = order.index(wanted)
            while j > i:
                path = [slots[t] for t in order]
                left, right = path[j-1:j+1]
                path[j-1:j+1] = [(left, right)]
                for term, coefficient in projection(tuple(path)).items():
                    out[term] += coefficient / factorial(n)
                order[j-1], order[j] = order[j], order[j-1]
                j -= 1
    return {t: c for t, c in out.items() if c}


def calculate():
    native = native26()
    tests, nonzero = Counter(), Counter()
    matrices = [fmpz_mat(a.tolist()) for a in native]
    for seed in (52071, 52072):
        rng = np.random.default_rng(seed)
        colours = [fmpz_mat(m.tolist()) for m in
                   (rng.integers(-2, 3, size=(4, 52)) @ native.reshape(52, -1)).reshape(4, 26, 26)]
        @cache
        def lie(node):
            if isinstance(node, int):
                return colours[node]
            x, y = map(lie, node)
            return x*y - y*x

        for n in range(1, 5):
            for slots in permutations(range(n)):
                terms = projection(slots)
                denominator = lcm(*(c.denominator for c in terms.values()))
                predicted = fmpz_mat(26, 26)
                for (value, pair), coeff in terms.items():
                    trace = trace_pair(lie(pair[0]), lie(pair[1])) if pair else 1
                    predicted += int(coeff*denominator)*trace*lie(value)
                actual = colours[slots[0]]
                for i in slots[1:]:
                    actual = actual * colours[i]
                difference = denominator * actual - predicted
                assert all(trace_pair(t, difference) == 0 for t in matrices)
                tests[n] += len(matrices)
                nonzero[n] += int(not predicted.is_zero())
        lie.cache_clear()
    shape = (2, 2, 2)
    values = [character(shape, cycles) for cycles in ((1,)*6, (1,)*4+(2,), (1,)*3+(3,))]
    assert values[0] + 3*values[1] + 2*values[2] == 0
    # Degree five first permits an adjoint covariant independent of X on the Cartan.
    x = sum((v*matrices[i] for i, v in enumerate((2, 3, 5, 7))), fmpz_mat(26, 26))
    cartan_columns = fmpz_mat([[trace_pair(t, x), trace_pair(t, x**5)] for t in matrices[:4]])
    assert cartan_columns.rank() == 2
    return dict(native_metric_convention="B=Tr26/6",inverse_metric_edge_factor=6,
                projection_test_components=dict(tests), nonzero_ordered_probes=dict(nonzero),
                quartic_polarization="1/36",
                flavour_S6_character_identity_transposition_three_cycle=values,
                invariant_dimension_for_one_symmetric_triple=0,
                cartan_probe=[2,3,5,7],
                linear_and_quintic_cartan_pairings=[[int(c) for c in r] for r in cartan_columns.tolist()])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(), indent=2))
