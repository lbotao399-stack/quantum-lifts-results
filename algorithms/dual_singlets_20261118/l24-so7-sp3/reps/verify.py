"""Exact rational checks in the universal signed dihedral trace algebra."""
import gzip
import json
from fractions import Fraction
from pathlib import Path

from trace_q0 import add, canonical, count, cyc, move, op

ROOT = Path(__file__).resolve().parent


def as_tuple(x):
    return tuple(map(as_tuple, x)) if isinstance(x, list) else x


def read(filename, algebra, weight, length, nterms=None):
    with gzip.open(ROOT / filename, "rt", encoding="utf-8") as f:
        raw = json.load(f)
    assert raw["algebra"] == algebra
    assert tuple(raw["weight"]) == weight and raw["Y"] == length
    terms = {as_tuple(w): Fraction(c) for w, c in raw["terms"]}
    assert len(terms) == len(raw["terms"])
    if nterms is not None:
        assert len(terms) == nterms
    for monomial, coefficient in terms.items():
        letters = sum(monomial, ())
        assert coefficient and len(letters) == length and count(letters) == weight
        assert canonical(monomial, op) == (1, monomial)
        assert all(cyc(word) == (word, 1) for word in monomial)
    return terms


def act(terms, k):
    result = {}
    for w, c in terms.items():
        add(result, move(w, k), c)
    return result


def main():
    P = read("sp3-p0.json.gz", "Sp(3)", (5, 2, 2, 1, 1), 9, 11536)
    R = read("sp3-p0-defect.json.gz", "Sp(3)", (5, 2, 2, 1, 1), 10, 16620)
    D = read("so7-cartan-null.json.gz", "SO(7)", (6, 3, 3, 0, 0), 10)
    for terms in (P, R, D):
        for k in range(3):
            assert not act(terms, k), k
    assert act(P, 3) == R
    assert not act(R, 3)
    assert not act(D, 3)
    print(f"PASS: compact highest weights; Q0 P = R; Q0 R = Q0 D = 0")
    print(f"Terms: Sp(3) P={len(P)}, defect R={len(R)}, SO(7) D={len(D)}")


if __name__ == "__main__":
    main()
