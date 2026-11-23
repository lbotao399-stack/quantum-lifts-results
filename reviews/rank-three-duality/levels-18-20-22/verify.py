"""Return the raising words of the SO(7) level-18 pair to their lowest weight."""
from collections import Counter, defaultdict
from functools import cache
from hashlib import sha1
from itertools import combinations, permutations, product
import json
from pathlib import Path

from jet_algebra import E, P, bracket, raising, sign, transpose, verify_brackets

HERE = Path(__file__).resolve().parent
WEIGHTS = ((1, 1, 2, 2, 2), (0, 0, 3, 3, 3))
CENSUS_BLOBS = {18: '8eacd8933cfbd6bae4fbd41b20726f41003a84f3',
                20: '0cc3ae93a62e5922132a6ea5a71ebf3ea1b916e2',
                22: '6b398d8d421358dd6f37015272000757202bd5de'}


def words(degree):
    """Poincare-Birkhoff-Witt basis in translations and linear odd raisings."""
    shears = [E(d, f) for d in range(2) for f in range(2, 5)]
    odd = [P(f) for f in range(2, 5)]
    for n in range(min(6, degree) + 1):
        for m in range(min(3, (degree - n) // 2) + 1):
            rest = degree - n - 2*m
            if rest % 3:
                continue
            for p in range(rest // 3 + 1):
                even = (P(0),) * p + (P(1),) * (rest // 3 - p)
                for theta, linear in product(combinations(odd, m), combinations(shears, n)):
                    yield even + theta + linear


def shift(word):
    w = [0] * 5
    for op in word:
        w[op[1]] += 1
        if op[0] == 'E':
            w[op[2]] -= 1
    return tuple(w)


def coefficient(primary):
    @cache
    def returned(word):
        if not word:
            return 1
        last = word[-1]
        if not raising(last):
            if last[0] == 'E' and last[1] == last[2]:
                return primary[last[1]] * returned(word[:-1])
            return 0
        for k in reversed(range(len(word)-1)):
            a, b = word[k:k+2]
            if not raising(a):
                start, tail = word[:k], word[k+2:]
                return (sign(a, b) * returned(start + (b, a) + tail)
                        + sum(v * returned(start + (op,) + tail) for v, op in bracket(a, b)))
        return 0
    return returned


def determinant(matrix):
    """Bareiss elimination, including the pivot sign."""
    a = [row.copy() for row in matrix]
    last, phase = 1, 1
    for j in range(len(a) - 1):
        swap = next((k for k in range(j, len(a)) if a[k][j]), None)
        if swap is None:
            return 0
        if swap != j:
            a[j], a[swap] = a[swap], a[j]
            phase = -phase
        for k in range(j+1, len(a)):
            for l in range(j+1, len(a)):
                n = a[j][j]*a[k][l] - a[k][j]*a[j][l]
                a[k][l], rem = divmod(n, last)
                assert rem == 0
        last = a[j][j]
    return phase * a[-1][-1]


def charges(weight):
    a, b, i, j, k = weight
    return (2*a+i+j+k-8, 2*b+i+j+k-8, a+b+2*i, a+b+2*j, a+b+2*k)


def highest(components):
    """SU(3) Weyl numerator at fixed Lorentz weights, in doubled charges."""
    result = {}
    for sector in components:
        spin, weight = sector[:2], sector[2:]
        if list(weight) != sorted(weight, reverse=True):
            continue
        count = 0
        for perm in permutations(range(3)):
            phase = (-1)**sum(x>y for k, x in enumerate(perm) for y in perm[k+1:])
            probe = tuple(x + 2*(perm[i]-i) for i, x in enumerate(weight))
            count += phase * components[spin+probe]
        if count:
            assert count > 0
            result[sector] = count
    return result


def run(degree):
    level = 18 + degree
    blocks = defaultdict(list)
    for w in words(degree):
        blocks[shift(w)].append(w)
    pairings = [coefficient(p) for p in WEIGHTS]
    records = []
    components = [Counter(), Counter()]
    for delta, basis in sorted(blocks.items()):
        dual = [tuple(transpose(op) for op in reversed(w)) for w in basis]
        matrices = [[[row(r+w) for w in basis] for r in dual] for row in pairings]
        assert matrices[0] == matrices[1]
        matrix = matrices[0]
        det = determinant(matrix)
        assert det != 0, (level, delta, matrix)
        for p, bank in zip(WEIGHTS, components):
            bank[charges(tuple(u+v for u, v in zip(p, delta)))] += len(basis)
        records.append(dict(shift=delta, words=basis, matrix=matrix, determinant=det))
    for row in pairings:
        row.cache_clear()
    counts = [highest(c) for c in components]
    joins = Counter(counts[0]) + Counter(counts[1])
    census = (HERE/f'census_j{level}.json').read_bytes()
    assert sha1(f'blob {len(census)}\0'.encode()+census).hexdigest() == CENSUS_BLOBS[level]
    source = json.loads(census)
    observed = {tuple(v['sec']): v['SO_odd']-v['Sp'] for v in source['differences']}
    assert all(v['Y'] == 8 for v in source['differences'])
    assert joins == observed
    return dict(level=level, component_rank=sum(len(v) for v in blocks.values()),
                multiplet_rank=[sum(c.values()) for c in counts], largest_block=max(map(len, blocks.values())),
                positive_charge_differences=[[w, joins[w]] for w in sorted(joins)], blocks=records)


if __name__ == '__main__':
    assert __debug__
    result = dict(single_jet_brackets=verify_brackets(), levels=[run(k) for k in (0, 2, 4)])
    print(json.dumps(result, indent=2))
