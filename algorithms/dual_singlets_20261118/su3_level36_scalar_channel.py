"""Resolve the two C-even level-36 descendants by their quartic scalar restriction."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path

from su3_bmn_descendants import (LEVEL34_SHA, PREVIOUS, WORDS, highest, load_words,
                                 polynomial, read_operator, restrict_block, witness)
from match_threshold import canonical, paper_core
from verify import add_scaled


def cochain_dimension(sector, words):
    # In the adjoint of SU(2), the weights in units of the positive root are -1,0,1.
    # Subtracting the weight-one multiplicity removes every nontrivial integer spin.
    dimension = 0
    for content in words.contents_of_sector(sector):
        character = {0: 1}
        for letter, count in enumerate(content):
            combinations_for_letter = combinations if words.is_ferm(letter) else combinations_with_replacement
            weights = Counter(sum(parts) for parts in combinations_for_letter((-1, 0, 1), count))
            updated = Counter()
            for old, multiplicity in character.items():
                for new, number in weights.items():
                    updated[old + new] += multiplicity * number
            character = updated
        dimension += character[0] - character[1]
    return dimension


def image_rank(sector, words, prime):
    pivots = {}
    for word in words.words_of_sector(sector, 3):
        row = {m: c % prime for m,c in polynomial(words.Q_multi(word)).items() if c % prime}
        while row:
            key = min(row)
            coefficient = row[key]
            if key not in pivots:
                inverse = pow(coefficient, -1, prime)
                pivots[key] = {m: c * inverse % prime for m,c in row.items()}
                break
            for m,c in pivots[key].items():
                reduced = (row.get(m, 0) - coefficient * c) % prime
                if reduced:
                    row[m] = reduced
                else:
                    row.pop(m, None)
    return len(pivots)


def classical_block(sector, words, prime):
    before = (sector[0]+1, *sector[1:])
    dimension = cochain_dimension(sector, words)
    incoming = image_rank(before, words, prime)
    outgoing = image_rank(sector, words, prime)
    return dict(sector=sector, cochains=dimension, incoming=incoming, outgoing=outgoing,
                rational_upper=dimension - incoming - outgoing)


def published_core(words):
    # arXiv:2412.08695v2 (2.6) uses the traceless bilinear and a normalized
    # symmetrization of the three fermion flavours. Each v_gen carries a factor 3.
    result = {}
    for p, q, r in permutations(range(3)):
        epsilon = (-1)**(int(p > q) + int(p > r) + int(q > r))
        for m, n in product(range(3), repeat=2):
            triple = words.sym_trace_expr((3 + m, 3 + n, 3 + p))
            pairs = words.mul_expr(words.v_gen((m,), q), words.v_gen((n,), r))
            add_scaled(result, words.mul_expr(pairs, triple), Fraction(epsilon, 3*3*6))
    return canonical(result, words)


def calculate(words, level34):
    prime = 1000003
    singlet_weight = classical_block((5, 4, 4, 4), words, prime)
    root_weight = classical_block((5, 5, 4, 3), words, prime)
    assert (singlet_weight["rational_upper"], root_weight["rational_upper"]) == (1, 0)
    core = published_core(words)
    assert polynomial(core) == polynomial(paper_core(words))
    core_record = witness(core, (5, 4, 4, 4), words)
    assert core_record["rational_class_nonzero"]
    # Apply Q^3_+ first, then Q^2_+, using the order of the cyclic epsilon contraction.
    double = words.Qplus_expr(words.Qplus_expr(core, 2), 1)
    double_record = witness(double, (7, 4, 5, 5), words)
    assert double_record["rational_class_nonzero"]
    raised = highest(words.Qplus_expr(level34, 2), words)
    pure = restrict_block(raised, (4, 0, 0), words)
    pure_record = witness(pure, (7, 4, 5, 5), words)
    mixed = restrict_block(raised, (3, 0, 0), words, 5)
    mixed_record = witness(mixed, (6, 4, 4, 5), words)
    assert not pure_record["rational_class_nonzero"] and mixed_record["rational_class_nonzero"]
    # The singlet coefficient of the level-32 parent is the input from arXiv:2412.08695v2.
    coefficient = -1020600
    pairing = coefficient * Fraction(double_record["pairing"])
    determinant = pairing * Fraction(mixed_record["pairing"])
    assert determinant
    return dict(prime=prime, weight=[8, 5, 5], length=11, colour_parity="even",
                active_level24=[singlet_weight, root_weight], core=core_record,
                twice_raised_core=double_record, level34_pure_four_scalar=pure_record,
                level34_mixed=mixed_record, published_level32_core_coefficient=coefficient,
                level32_pure_pairing=str(pairing), determinant=str(determinant))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--words", type=Path, default=WORDS)
    parser.add_argument("--level34", type=Path, default=PREVIOUS/"O34_C_even.json.gz")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    record = calculate(load_words(args.words), read_operator(args.level34, LEVEL34_SHA))
    args.out.write_text(json.dumps(record, indent=2) + "\n")
