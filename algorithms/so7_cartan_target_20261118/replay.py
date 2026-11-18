"""Restrict the saved SO(7) traces to three independent Cartan planes over Q.

Each letter on a plane is its own commuting or exterior variable times the
matrix ((0,1),(-1,0)). The seventh color has zero weight. A nonempty trace
therefore contributes 2*(-1)**(length/2) on each plane of even length.
The second evaluation uses diagonal weights in opposite pairs and multiplies
one trace polynomial at a time using exterior bit masks.
"""
import gzip
import hashlib
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import lcm, prod
from pathlib import Path


def trace_expansion(expression, odd_letters):
    denominator = lcm(*(row[1] for row in expression))
    total = Counter()
    contributing_words = 0
    for numerator, divisor, word in expression:
        if any(len(trace) % 2 for trace in word):
            continue
        trace_factor = prod(2 * (-1) ** (len(trace) // 2) for trace in word)
        one_word = Counter()
        for planes in product(range(3), repeat=len(word)):
            even, odd = [], []
            for trace, plane in zip(word, planes):
                for letter in trace:
                    (odd if letter in odd_letters else even).append(3 * letter + plane)
            if len(odd) != len(set(odd)):
                continue
            inversions = sum(x > y for i, x in enumerate(odd) for y in odd[i + 1:])
            key = tuple(sorted(even)), tuple(sorted(odd))
            one_word[key] += (-1) ** inversions * trace_factor
        nonzero = {k: v for k, v in one_word.items() if v}
        contributing_words += bool(nonzero)
        for key, value in nonzero.items():
            total[key] += numerator * (denominator // divisor) * value
    polynomial = {k: Fraction(v, denominator) for k, v in total.items() if v}
    return polynomial, {
        "input_words": len(expression),
        "coefficient_denominator": denominator,
        "nonzero_words_before_combination": contributing_words,
        "distinct_monomials_before_combination": len(total),
        "nonzero_monomials": len(polynomial),
    }


def multiply(left, right):
    result = Counter()
    for (a, mask), x in left.items():
        for (b, new), y in right.items():
            if mask & new:
                continue
            bits, sign = new, 1
            while bits:
                bit = bits & -bits
                if (mask >> bit.bit_length()).bit_count() % 2:
                    sign = -sign
                bits ^= bit
            result[tuple(sorted(a + b)), mask | new] += sign * x * y
    return {k: v for k, v in result.items() if v}


def diagonal_evaluator(odd_letters):
    unit = {((), 0): 1}

    @lru_cache(None)
    def trace_polynomial(trace):
        result = Counter()
        for plane, weight in product(range(3), (-1, 1)):
            term = unit
            for letter in trace:
                index = 3 * letter + plane
                key = ((), 1 << index) if letter in odd_letters else ((index,), 0)
                term = multiply(term, {key: weight})
            result.update(term)
        return {k: v for k, v in result.items() if v}

    def evaluate(expression):
        result = Counter()
        for numerator, divisor, word in expression:
            term = unit
            for trace in word:
                term = multiply(term, trace_polynomial(tuple(trace)))
                if not term:
                    break
            for key, value in term.items():
                result[key] += Fraction(numerator, divisor) * value
        return {k: v for k, v in result.items() if v}

    return evaluate


def replay():
    blob = Path(__file__).with_name("representatives.json.gz").read_bytes()
    assert hashlib.sha256(blob).hexdigest() == "eda182846ea9ac1f2ed989519f5b4a234e35ccf48f0ed6f120ce899b0ac37e20"
    data = json.loads(gzip.decompress(blob))
    odd = {i for i, letter in enumerate(data["letters"]) if len(letter[2]) % 2 == 0}
    diagonal = diagonal_evaluator(odd)
    evaluated, rows = {}, {}
    for name in ("S", "F", "SF"):
        expr = data[name]
        lengths = {sum(map(len, word)) for _, _, word in expr}
        assert len(lengths) == 1
        length = lengths.pop()
        assert length % 2 == 0
        direct, rows[name] = trace_expansion(expr, odd)
        evaluated[name] = diagonal(expr)
        # Eigenvalues of a skew block are +/- i times its coefficient.
        common = {(even, sum(1 << i for i in exterior)): value
                  for (even, exterior), value in direct.items()}
        assert common == {k: (-1) ** (length // 2) * v for k, v in evaluated[name].items()}
    assert len(evaluated["F"]) == 153
    assert rows["SF"]["distinct_monomials_before_combination"] == 21
    assert multiply(evaluated["S"], evaluated["F"]) == evaluated["SF"] == {}
    square = multiply(evaluated["S"], evaluated["S"])
    cube = multiply(square, evaluated["S"])
    assert len(cube) == 1 and not multiply(cube, evaluated["S"])
    return {"group": "SO(7)", "level": 30, "field": "Q",
            "input_certificate_sha256": data["source_certificate_sha256"],
            "polynomials": rows, "diagonal_weight_crosscheck": True,
            "factorized_product_crosscheck": True, "S_cubed_nonzero_on_Cartan": True,
            "S_fourth_zero_on_Cartan": True}


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2))
