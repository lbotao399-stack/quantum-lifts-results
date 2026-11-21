"""Recover rational chart coordinates at a common denominator using exact equation rows."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd

from su3_subregular_differential import ODD, qmonomial, variation


@lru_cache(maxsize=1)
def reverse_variations():
    reverse = {}
    for index in range(len(ODD)):
        if index >= 64 and index % 16 == 0:
            continue
        for word in variation(index):
            reverse.setdefault(word, []).append(index)
    return reverse


def support_row(target, residues, image):
    parents = set()
    for length in (2, 3):
        for positions in combinations(range(len(target)), length):
            terms = reverse_variations().get(bytes(target[i] for i in positions), ())
            if not terms:
                continue
            rest = bytes(v for i, v in enumerate(target) if i not in positions)
            for k in terms:
                word = bytes(sorted((*rest, k)))
                if word in residues:
                    parents.add(word)
    return [(residues[w], n) for w in sorted(parents) if (n := image(w).get(target))]


def recover(residues, modulus, unit, rhs, degree, image_of, *, limit=65536):
    """Try bounded common denominators. Only an exact multiplication can return a solution."""
    if image_of is None:
        def image(w):
            return {v: n for v, n in qmonomial(w).items() if len(v) == degree}
    else:
        image = image_of

    integer_rhs = {}
    for word, coefficient in rhs.items():
        value = coefficient / unit
        if value.denominator != 1:
            raise ValueError("The exported right hand side has a nonintegral coordinate")
        integer_rhs[word] = value.numerator

    # Every failed exact replay supplies further rows that the next denominator must satisfy.
    # The centred numerator may wrap at a small modulus. Failure only requests another prime.
    half = modulus // 2
    tests = []
    trials = 0
    for denominator in range(1, limit + 1):
        if gcd(denominator, modulus) != 1:
            continue
        if any(sum(c * ((r * denominator + half) % modulus - half) for r, c in row) != denominator * b
               for row, b in tests):
            continue
        trials += 1
        numerators = {w: (r * denominator + half) % modulus - half for w, r in residues.items() if r}
        residual = {w: denominator * c for w, c in integer_rhs.items()}
        for w, value in numerators.items():
            for t, c in image(w).items():
                remaining = residual.get(t, 0) - value*c
                if remaining:
                    residual[t] = remaining
                else:
                    residual.pop(t, None)
        if not residual:
            return {w: unit * Fraction(n, denominator) for w, n in numerators.items()}, denominator, trials

        # Several well-separated failing rows distinguish multiples of an earlier denominator.
        step = max(1, len(residual) // 24)
        selected = sorted(residual)[::step][:24]
        for target in selected:
            row = support_row(target, residues, image)
            b = integer_rhs.get(target, 0)
            error = denominator*b - sum(c*((r*denominator+half) % modulus-half) for r,c in row)
            if error != residual[target]:
                raise ArithmeticError("Backward and forward sparse coefficients disagree")
            tests.append((row, b))
        if trials == 12:
            break
    return None
