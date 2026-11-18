"""Identify the mixed block component modulo the full SU(2) boundary space."""

import gzip
import hashlib
import json
from collections import Counter
from fractions import Fraction
from itertools import permutations, product

from verify import (HERE, add_scaled, boundary_basis, even_part, load_source, load_words,
                    pairing, polynomial, reduce_row, restrict_block)


def paper_core(words):
    """Kim, Lee, Lee and Song, equation (3.6), with the three traces in their printed order."""
    source = Counter()
    for i, j, k in permutations(range(3)):
        sign = (-1) ** (int(i > j) + int(i > k) + int(j > k))
        for m, n in product(range(3), repeat=2):
            factor = sign
            factors = []
            for trace in ((m, 3+i), (n, 3+j), (3+m, 3+n, 3+k)):
                cyclic, representative = words.canon_word(trace)
                if not cyclic:
                    break
                factors.append(representative)
                factor *= cyclic
            else:
                exchange, word = words.canon_multi(factors)
                if exchange:
                    source[word] += factor * exchange
    return {word: value for word, value in source.items() if value}


def canonical(expression, words):
    result = Counter()
    for traces, coefficient in expression.items():
        factors = []
        for trace in traces:
            cyclic, representative = words.canon_word(trace)
            if not cyclic:
                break
            coefficient *= cyclic
            factors.append(representative)
        else:
            exchange, word = words.canon_multi(factors)
            if exchange:
                result[word] += coefficient * exchange
    return {w:c for w,c in result.items() if c}


def compare():
    words = load_words(HERE.parent / "su3_independent_prime_20261118"
                       / "runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor/fortuity/bmn/words.py")
    original = load_source(HERE / "published_O34_input.json.gz", words)
    source = even_part(original, words)
    raw = (HERE / "odd_parity_primitive.json.gz").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "3af11db9cc1d8f47ce29239741ac3fe99833838871de8e88a84fec17ccf93c1c"
    parity_primitive = {}
    for traces, c in json.loads(gzip.decompress(raw)):
        word = tuple(tuple(6 if label == 7 else label for label in t) for t in traces)
        parity_primitive[word] = Fraction(c)
    odd = canonical(original, words)
    add_scaled(odd, source, -1)
    odd_terms = len(odd)
    add_scaled(odd, words.Q_expr(canonical(parity_primitive, words)), Fraction(-1,2))
    assert not odd
    image = polynomial(restrict_block(source, (3, 0, 0), words, 5))
    mixed = json.loads((HERE / "certificate.json").read_text())["channels"][2]
    ell = {(e, o): Fraction(c) for e, o, c in mixed["witness"]["monomial_coefficients"]}
    core = paper_core(words)
    assert not polynomial(words.Q_expr(core))
    core_value = polynomial(core)
    core_pairing = pairing(ell, core_value)
    source_pairing = pairing(ell, image)
    assert (source_pairing, core_pairing) == (120294720, -20)
    multiple = Fraction(source_pairing, core_pairing)
    difference = image.copy()
    add_scaled(difference, core_value, -multiple)
    predecessors = words.words_of_sector((6, 4, 4, 4), 3)
    basis = boundary_basis(predecessors, words)
    remainder, used = reduce_row(difference, basis)
    primitive = {predecessors[i]: c for i, c in used.items()}
    assert not remainder and polynomial(words.Q_expr(primitive)) == difference
    return dict(source="Kim, Lee, Lee and Song, arXiv:2610.20340v2, equation (3.6)",
                parity_predecessor_terms=len(parity_primitive), parity_predecessor_coefficient="1/2",
                odd_source_terms=odd_terms, formal_parity_residual_terms=len(odd),
                core_trace_terms=len(core), core_monomials=len(core_value),
                mixed_pairing=str(source_pairing), core_pairing=str(core_pairing),
                mixed_class_over_core=str(multiple), difference_monomials=len(difference),
                primitive_trace_terms=len(primitive), residual_monomials=len(remainder),
                primitive=[[w, str(c)] for w, c in sorted(primitive.items())])


if __name__ == "__main__":
    record = compare()
    assert json.loads(json.dumps(record)) == json.loads((HERE / "threshold_identity.json").read_text())
    print(json.dumps(record, indent=2))
