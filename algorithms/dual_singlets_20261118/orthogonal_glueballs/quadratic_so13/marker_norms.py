"""Check the scalar determinant factors, the strict charges and both families of first pairings."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, prod
import json
from auxiliary import HERE
from verify import target_coordinates
from raising import raw_word
from so8_level24_descendant import check_letter, substitute, substitution
from q1d_exact import pair_terms
from rank_four_lowered_replay import terminal


def norm(poly):
    return sum(
        c * c * prod(factorial(k) for k in Counter(w).values()) for w, c in poly.items()
    )


def times(left, right):
    out = Counter()
    for w, c in left.items():
        for v, b in right.items():
            out[tuple(sorted(w + v))] += c * b
    return out


def run():
    det = {
        tuple(3 * i + p[i] for i in range(3)): (-1)
        ** sum(a > b for i, a in enumerate(p) for b in p[i + 1 :])
        for p in permutations(range(3))
    }
    minor = {(1, 5): 1, (2, 4): -1}
    assert [norm(det), norm(times(det, minor)), norm(times(det, det))] == [6, 24, 144]
    jets = {
        n: [
            (i, j, t)
            for i in range(n + 1)
            for j in range(n - i + 1)
            for t in combinations(range(3), n - i - j)
        ]
        for n in (1, 2)
    }
    pairs = (
        list(product(jets[1], jets[1]))
        + list(product(jets[1], jets[2]))
        + list(product(jets[2], jets[1]))
    )
    for i in range(3):
        for letter in jets[1] + jets[2]:
            check_letter(letter, 1, i, True)
        for x, y in pairs:
            for z, k in ((x, 0), (y, 1)):
                c, v = substitution(z, 1, i, True)
                if c:
                    assert not (pair_terms(v, y) if k == 0 else pair_terms(x, v))
            for _, a, b in pair_terms(x, y):
                assert (
                    not substitution(a, 1, i, True)[0]
                    and not substitution(b, 1, i, True)[0]
                )
    for m, file, factor in ((2, "triplet", 4), (3, "quartet", 24)):
        for n, p, v in json.loads((HERE / (file + ".json")).read_text()):
            k = p + m
            power_of_two = Fraction(2) ** (k - 5)
            assert (
                Fraction(v) ** 2
                == (factor * 3 * power_of_two * (2 * n - 3 * k - 4)) ** 2
            )
            basis, _ = target_coordinates(
                (p + m, p) + tuple(int(i < m) for i in range(3))
            )
            assert len(basis) == 1
            poly = raw_word(basis[0])
            for i in range(m):
                poly = substitute(poly, 1, i, True)
            assert abs(sum(c * terminal(w)[0] for w, c in poly.items())) == 1
    print(
        "Marker norms 6, 24, 144. All 32 signed entries pass the squared formula and the strict target substitution."
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    run()
