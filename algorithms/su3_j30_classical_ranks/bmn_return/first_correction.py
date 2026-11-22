"""Reconstruct the first correction to the nine-letter BMN return from two-trace equations."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import product
from flint import fmpq_mat
from fortuity.full.q1 import Q1_expr
from fortuity.full.words import canon_multi, word_parity


def multiply(*rows):
    result = {(): 1}
    for row in rows:
        expanded = Counter()
        for left, a in result.items():
            for right, b in row.items():
                sign, word = canon_multi(left + right)
                if sign:
                    expanded[word] += sign * a * b
        result = {word: c for word, c in expanded.items() if c}
    return result


def parents(word):
    for i, trace in enumerate(word):
        if len(trace) < 3:
            continue
        for j in range(len(trace)):
            a, b, *rest = trace[j:] + trace[:j]
            if set(a[2]).intersection(b[2]):
                continue
            merged = (a[0]+b[0], a[1]+b[1], tuple(sorted(a[2]+b[2])))
            sign, source = product.su.canonical_multi(word[:i] + ((merged, *rest),) + word[i+1:])
            if sign:
                yield source


def primitive(rhs):
    """First try parents of the support, then propagate incidence only if needed."""
    if not rhs:
        return {}, 0
    frontier = set(rhs)
    seen = set(rhs)
    sources = {}
    while frontier:
        new = {p for w in frontier for p in parents(w)} - sources.keys()
        if not new:
            break
        sources.update((p,product.su.differential({p:1})) for p in sorted(new))
        words = sorted(rhs.keys() | {w for col in sources.values() for w in col})
        idx = {w:i for i,w in enumerate(words)}
        columns = sorted(sources)
        matrix = fmpq_mat(len(words), len(columns)+1)
        for j,p in enumerate(columns):
            for w,c in sources[p].items():
                matrix[idx[w],j] = c
        for w,c in rhs.items():
            matrix[idx[w],len(columns)] = str(-c)
        reduced, rank = matrix.rref()
        answer = {}
        for i in range(rank):
            j = next(j for j in range(len(columns)+1) if reduced[i,j])
            if j == len(columns):
                break
            value = Fraction(str(reduced[i,len(columns)]))
            if value:
                answer[columns[j]] = value
        else:
            assert product.su.differential(answer) == {w:-c for w,c in rhs.items()}
            return answer, len(columns)
        frontier = set(words)-seen
        seen.update(frontier)
    raise ValueError("The pair has a nonzero residual in the free trace space")


def encoded(poly):
    return [[w,str(c)] for w,c in sorted(poly.items())]


def construct():
    factors = [product.trace_factor(k,unitary=True) for k in product.FACTORS]
    parities = []
    for factor in factors:
        parity, = {word_parity(tr) for (tr,) in factor}
        parities.append(parity)
        assert not product.su.differential(factor) and not Q1_expr(factor,3)
    raw = Counter()
    records = []
    for i,j in combinations(range(4),2):
        pair = multiply(factors[i],factors[j])
        rhs = Q1_expr(pair,3)
        correction, parent_count = primitive(rhs)
        rest = [k for k in range(4) if k not in (i,j)]
        inversions = sum(parities[k]*parities[l] for k in (i,j) for l in rest if l<k)
        for w,c in multiply(correction,*(factors[k] for k in rest)).items():
            raw[w] += (-1)**inversions*c
        records.append(dict(factors=[i+1,j+1],image_terms=len(rhs),parents_tested=parent_count,
                            primitive=encoded(correction)))
    raw = {w:c for w,c in raw.items() if c}
    unprojected = multiply(*factors)
    assert product.su.differential(raw) == {w:-c for w,c in Q1_expr(unprojected,3).items()}
    source = product.bmn_product()
    correction = product.su.projection.project(raw, product.sector((1,1,4,3,3),8))
    correction = {w:Fraction(c,2) for w,c in product.su.parity_project(correction,-1).items()}
    rhs = Q1_expr(source,3)
    assert product.su.differential(correction) == {w:-c for w,c in rhs.items()}
    assert product.charge(correction) == ((1,1,4,3,3),8)
    assert product.su.conjugate(correction) == {w:-c for w,c in correction.items()}
    for root in ((0,1),(1,2),("spin",1)):
        assert not product.su.projection.action(correction,*root)
    numbers,denominator = product.integral(correction)
    receipt = dict(source_terms=len(source),single_factor_images=[0]*4,pairs=records,
                raw_correction_terms=len(raw),projected_correction_terms=len(correction),
                projected_image_terms=len(rhs),residual_terms=0,correction_denominator=denominator,
                correction_numerator_sha256=product.digest(numbers))
    return source,numbers,receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    result = json.dumps(construct()[2],indent=2)+"\n"
    if args.check:
        assert args.check.read_text() == result
    else:
        print(result,end="")
