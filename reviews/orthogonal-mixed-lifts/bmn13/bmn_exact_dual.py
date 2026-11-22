"""Rebuild the rational covector with a symmetric-group Casimir projector."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial, gcd, lcm, prod
from pathlib import Path
import sys

from bmn_native_pairs import LETTERS, PRIMES, decode, load_duals


def partitions(total, bound, rows):
    if not total:
        yield ()
    elif rows:
        for i in range(min(total, bound), 0, -1):
            for tail in partitions(total-i, i, rows-1):
                yield (i,)+tail


def projector():
    highest = (4, 2, 2)
    weight = (0,)*4+(1,)*2+(2,)*2
    shapes = [shape for shape in partitions(8, 8, 3)
              if all(sum(shape[:i]) >= sum(highest[:i]) for i in (1, 2))]
    content = lambda shape: sum(j-i for i, n in enumerate(shape) for j in range(n))
    eigenvalue = content(highest)
    assert all(content(s) > eigenvalue for s in shapes if s != highest)
    others = sorted({content(s) for s in shapes if s != highest})
    result, divisor = {weight: 1}, 1
    for eigen in others:
        row = Counter({w: -eigen*v for w, v in result.items()})
        for w, v in result.items():
            for i, j in combinations(range(8), 2):
                u = list(w)
                u[i], u[j] = u[j], u[i]
                row[tuple(u)] += v
        result = {w: v for w, v in row.items() if v}
        divisor *= eigenvalue-eigen
    common = gcd(abs(divisor), *map(abs, result.values()))
    common *= 1 if divisor > 0 else -1
    return {w: v//common for w, v in result.items()}, divisor//common, others


@lru_cache(32768)
def stabilizer(trace):
    orbit = {t[i:]+t[:i] for t in (trace, trace[::-1]) for i in range(len(t))}
    return 2*len(trace)//len(orbit)


def metric(w):
    return prod(map(stabilizer, w))*prod(factorial(v) for v in Counter(w).values())


def run(a):
    sys.path.insert(0, str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from mixed_incoming import so_canon
    from so8_level24_descendant import substitute
    signed = {}
    for line in a.seed.read_text().splitlines():
        value, raw = line.split()
        s, w = so_canon(decode(raw))
        assert s and w not in signed
        signed[w] = Fraction(s*int(value), metric(w))
    den = lcm(*(v.denominator for v in signed.values()))
    coefficients, divisor, spectrum = projector()
    projected = Counter()
    for w, v in signed.items():
        positions = sorted((l[2][0], i, j) for i, tr in enumerate(w)
                           for j, l in enumerate(tr) if l[2])
        assert tuple(f for f, _, _ in positions) == (0,)*4+(1, 1, 2, 2)
        tt = [list(t) for t in w]
        scale = int(v*den)
        for assignment, k in coefficients.items():
            for flavour, (_, i, j) in zip(assignment, positions):
                tt[i][j] = LETTERS[flavour]
            sign, key = so_canon(tuple(map(tuple, tt)))
            if sign:
                projected[key] += sign*scale*k
    projected = {w: v for w, v in projected.items() if v}
    counts = [len(projected)]
    for spin, flavour in reversed([(i, j) for i in range(2) for j in range(3)]):
        projected = substitute(projected, spin, flavour, False)
        counts.append(len(projected))
    rational = {w: Fraction(k*metric(w), den*divisor) for w, k in projected.items()}
    native = load_duals(a.modular)
    assert rational.keys() == native.keys()
    for w, v in rational.items():
        for p, expected in zip(PRIMES, native[w]):
            assert (v.numerator-expected*v.denominator) % p == 0, (w, v, expected, p)
    enc = {l: str(i) for i, l in enumerate(LETTERS)}
    with a.out.open('w') as f:
        for w, v in sorted(rational.items()):
            print(v, '/'.join(''.join(enc[l] for l in t) for t in w), file=f)
    print(dict(seed_den=den, content_spectrum=[2, *spectrum], projector_den=divisor,
               labelled_support=len(coefficients), stages=counts, dual_den=lcm(
                   *(v.denominator for v in rational.values()))), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('repo', 'seed', 'out'):
        p.add_argument('--'+k, type=Path, required=True)
    p.add_argument('--modular', type=Path, nargs=2, required=True)
    assert __debug__
    run(p.parse_args())
