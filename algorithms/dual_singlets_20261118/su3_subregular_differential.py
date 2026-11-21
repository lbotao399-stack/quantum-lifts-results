"""Check the subregular substitution against every signed classical Taylor coefficient."""
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path

from su3_subregular_image import FIELDS, ODD, THETA_DEGREE, coefficient, entries, multiply, ordered

RULES = (
    ((-2, "qfr"), (2, "qgu")), ((1, "fr"), (-1, "gu")),
    ((3, "au"), (-1, "pf"), (2, "qqf")), ((-3, "ar"), (1, "pg"), (-2, "qqg")),
    ((-2, "qfg"),), ((1, "fg"),), ((3, "af"),), ((-3, "ag"),),
)


@lru_cache(None)
def variation(index):
    field, jet = divmod(index, 16)
    sign = (-1)**THETA_DEGREE[jet % 4]
    values = Counter()
    for n, factors in RULES[field]:
        for term, c in coefficient(factors, jet).items():
            values[term] += sign*n*c
    return {m: n for m, n in values.items() if n}


def qmonomial(word):
    result = Counter()
    sign = 1
    for i, variable in enumerate(word):
        for new, c in variation(variable).items():
            left, prefix = ordered(word[:i], new)
            right, monomial = ordered(prefix, word[i+1:])
            result[monomial] += sign*left*right*c
        sign *= -1 if ODD[variable] else 1
    return {m: n for m, n in result.items() if n}


def qpoly(poly):
    result = Counter()
    for m, n in poly.items():
        for v, c in qmonomial(m).items():
            result[v] += n*c
    return {m: n for m, n in result.items() if n}


def code(letter):
    z, w, theta = letter
    return (2*z+w)*8 + sum(1 << i for i in theta)


def check():
    from su3_square_zero_image import ALPH
    from fortuity.full.letters import Q_letter

    indices = [i for i in range(128) if i < 64 or i % 16]
    assert all(not qpoly(variation(i)) for i in indices)
    rules = 0
    monomials = 0
    for letter in ALPH:
        rhs = [Counter() for _ in range(9)]
        for n, left, right in Q_letter(letter):
            rules += 1
            X, Y = entries(code(left)), entries(code(right))
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for m, c in multiply(X[3*i+k], Y[3*k+j], 4).items():
                            rhs[3*i+j][m] += n*c
        for poly, target in zip(entries(code(letter)), rhs):
            actual = qpoly(poly)
            assert actual == {m: n for m, n in target.items() if n}
            monomials += len(actual)
    return dict(chart_fields=FIELDS, component_generators=len(indices), square_zero_generators=len(indices),
                source_letters=len(ALPH), source_ordered_splits=rules, matrix_entries_checked=9*len(ALPH),
                nonzero_entry_monomials=monomials)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args=parser.parse_args()
    if not __debug__:
        parser.error("verification requires assertions")
    args.out.write_text(json.dumps(check(), indent=2)+"\n")
