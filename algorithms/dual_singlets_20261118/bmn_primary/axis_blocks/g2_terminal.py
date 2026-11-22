"""Exact six-first-jet colour character at two Taylor degrees (3,3) in G2."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import gcd

from flint import fmpz_mat

from transverse import colour_data, signed_orbit


@cache
def one_axis_cycles():
    weights, bracket = colour_data("G2")
    triples = defaultdict(list)
    for w in combinations(range(len(weights)), 3):
        triples[tuple(map(sum, zip(*(weights[i] for i in w))))].append(w)
    cycles = []
    for weight, words in sorted(triples.items()):
        matrix = fmpz_mat(14 * 14, len(words))
        for col, w in enumerate(words):
            for i, j in combinations(range(3), 2):
                for k, value in bracket[w[i], w[j]]:
                    matrix[14 * w[3 - i - j] + k, col] += (-1) ** (i + j - 1) * value
        kernel, count = matrix.nullspace()
        for col in range(count):
            terms = {w: int(kernel[i, col]) for i, w in enumerate(words) if kernel[i, col]}
            divisor = gcd(*terms.values())
            cycles.append((weight, [(w, value // divisor) for w, value in terms.items()]))
    assert len(cycles) == 182
    return cycles


@cache
def cross_image(i, j):
    cycles = one_axis_cycles()
    _, bracket = colour_data("G2")
    terms = Counter()
    for u, n in cycles[i][1]:
        for v, m in cycles[j][1]:
            for a, x in enumerate(u):
                for b, y in enumerate(v):
                    for z, c in bracket[x, y]:
                        left, right = u[:a] + u[a + 1:], v[:b] + v[b + 1:]
                        coefficient = (-1) ** (a + b) * n * m * c
                        terms[left, right, z] += coefficient
                        if i != j:
                            terms[right, left, z] -= coefficient
    assert all(value == -terms[right, left, z] for (left, right, z), value in terms.items())
    return {k: n for k, n in terms.items() if k[0] < k[1] and n}


@cache
def terminal_character():
    cycles = one_axis_cycles()
    by_weight = defaultdict(list)
    for i, j in combinations_with_replacement(range(len(cycles)), 2):
        weight = tuple(x + y for x, y in zip(cycles[i][0], cycles[j][0]))
        if min(weight) >= 0:
            by_weight[weight].append((i, j))

    result, checks = {}, []
    for weight, pairs in sorted(by_weight.items()):
        columns = [cross_image(i, j) for i, j in pairs]
        # The symmetric input gives an antisymmetric result on the two remaining pairs.
        # Choosing one orientation gives the same rank as keeping both.
        sites = sorted({k for c in columns for k in c})
        index = {w: j for j, w in enumerate(sites)}
        matrix = fmpz_mat(len(index), len(pairs))
        for col, values in enumerate(columns):
            for w, n in values.items():
                if w in index:
                    matrix[index[w], col] = n
        rank = matrix.rank()
        n = len(pairs) - rank
        checks.append([*weight, len(pairs), len(sites), rank, n])
        if n:
            result.update({w: n for w in signed_orbit(*weight, "G2")})
    return result, checks
