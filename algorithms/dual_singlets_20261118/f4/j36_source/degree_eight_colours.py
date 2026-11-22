"""Resolve degree-eight F4 invariant tensors into symmetric-group multiplicities."""
import argparse
from collections import Counter
from functools import cache
from itertools import permutations, product
import json
from random import Random

import numpy as np
from flint import nmod_mat

from mixed_colour_relations import (evaluate_words, identify, pivots, sample_colours,
                                    specialized_quartic_rank, trace_words)
from quartic_arm_complex import schur_multiplicity, sign
from b2_schur_complex import dimension, partitions
from f4_colour import basis, jordan_colours
from f4_singlet_characters import decomposition


@cache
def alternations(shape):
    columns = [[sum(shape[:row]) + col for row, size in enumerate(shape) if col < size]
               for col in range(shape[0])]
    choices = [tuple(permutations(range(len(col)))) for col in columns]
    result = []
    for order in product(*choices):
        labels, factor = [0] * sum(shape), 1
        for col, perm in zip(columns, order):
            factor *= sign(perm)
            for i, row in zip(col, perm):
                labels[i] = row
        result.append((tuple(labels), factor))
    return result


def projected_trace(word, shape):
    polynomial = Counter()
    for labels, factor in alternations(shape):
        key, orientation = identify(word, labels)
        polynomial[key] += factor * orientation
    return {k:v for k,v in polynomial.items() if v}


def quartic_pair(colours, jordan, inverse, prime):
    x, y = colours
    arms = []
    for k in range(5):
        word = [x]*(4-k) + [y]*k
        matrix = word[0]
        for colour in word[1:]:
            matrix = (matrix @ colour) % prime
        arms.append(np.einsum('tij,ji->t', jordan, matrix) % prime)
    def pairing(i, j):
        return int(((arms[i]*inverse % prime)*arms[j]).sum() % prime)
    return (2*pairing(0,4) - 8*pairing(1,3) + 6*pairing(2,2)) % prime


def block(shape, words, representation, jordan, inverse, prime):
    upper = schur_multiplicity(shape)
    if upper == 0:
        return dict(shape=shape, multiplicity=0, defining=0)
    expected = upper - int(shape == (4, 4))
    order = list(range(len(words)))
    Random(854321).shuffle(order)
    samples = [sample_colours(representation, len(shape), i, prime) for i in range(upper + 3)]
    selected = []
    rows = [[] for sample in samples]
    for start in range(0, len(order), 16):
        indices = order[start:start+16]
        polynomials = [projected_trace(words[j], shape) for j in indices]
        terms = sorted(set().union(*(p.keys() for p in polynomials)))
        if not terms:
            continue
        for row, sample in zip(rows, samples):
            values = dict(zip(terms, evaluate_words(terms, sample, prime)))
            row.extend(sum(v*values[t] for t,v in poly.items()) % prime for poly in polynomials)
        selected.extend(indices)
        columns = pivots(nmod_mat(rows, prime))
        if len(columns) == expected:
            break
        assert len(columns) < expected
    assert len(columns) == expected, (shape, len(columns), expected)
    matrix = [[row[j] for j in columns] for row in rows]
    if shape == (4, 4):
        for row, sample in zip(matrix, samples):
            row.append(quartic_pair(sample, jordan, inverse, prime))
    chosen = pivots(nmod_mat(matrix, prime).transpose())
    assert len(chosen) == upper
    determinant = int(nmod_mat([matrix[i] for i in chosen], prime).det())
    assert determinant
    return dict(shape=shape, multiplicity=upper, defining=expected, rows=chosen,
                trace_words=[words[selected[j]] for j in columns], determinant=determinant)


def run(prime):
    representation, _ = basis()
    jordan = jordan_colours()
    gram = np.einsum('iab,jba->ij', jordan, jordan)
    inverse = np.array([pow(int(v), -1, prime) for v in gram.diagonal()], dtype=np.int64)
    labelled = trace_words(8)
    records = [block(p[::-1], labelled, representation, jordan, inverse, prime)
               for p in partitions(8)]
    # Weight (4,4) has one vector in each two-row summand. Quartic
    # Lie identities bound the dimension of its entire defining-trace span.
    labels = (0,)*4 + (1,)*4
    words = sorted({t for w in labelled for t, s in [identify(w, labels)] if s})
    relations = specialized_quartic_rank(labelled, words, sorted(set(permutations(labels))), prime)
    assert len(words) - relations == sum(r['defining'] for r in records if len(r['shape']) <= 2)
    full = sum(dimension(r['shape'])*r['multiplicity'] for r in records)
    assert full == sum(n*n for n in decomposition('F4', (1,)*4).values())
    native = sum(dimension(r['shape'])*r['defining'] for r in records)
    return dict(prime=prime, full=full, defining=native, missing=full-native,
                specialized_trace_words=len(words), specialized_quartic_relations=relations,
                blocks=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, required=True, choices=(1000003, 1000033))
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    print(json.dumps(run(args.prime), separators=(',', ':')))
