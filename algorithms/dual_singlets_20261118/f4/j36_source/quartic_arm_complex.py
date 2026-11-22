"""Classical source and image of a contracted pair of quartic Albert covariants."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from flint import fmpz_mat, nmod_mat

from jordan_source import coefficients, colour_samples, evaluation_rows, columns
from mixed_colour_relations import pivots, representation as archive_matrices, trace_words
from verify import ARCHIVE, HERE
from b2_schur_complex import partitions, centralizer, character, multiplicity, transpose
from f4_colour import basis, jordan_colours
from f4_singlet_characters import ZERO, decomposition


@cache
def schur_multiplicity(shape):
    terms = (Fraction(character(shape, cycles) * decomposition('F4', cycles).get(ZERO, 0),
                      centralizer(cycles)) for cycles in partitions(sum(shape)))
    total = sum(terms)
    assert total.denominator == 1 and total >= 0
    return int(total)


def horizontal_strip(shape, size):
    for parts in partitions(sum(shape) + size):
        outer = parts[::-1]
        if len(outer) < len(shape):
            continue
        inner = shape + (0,) * (len(outer) - len(shape))
        if all(a >= b for a, b in zip(outer, inner)) and all(
                a >= b for a, b in zip(inner, outer[1:])):
            yield outer


def dimensions():
    # The arms are both closed when the double jet is in the tail.
    # Only the (4,4) part of Sym^2(Sym^4) fits inside the final diagram.
    # In the residual two Taylor directions the double jet transforms by det.
    def source_weight(scalars):
        result = Fraction()
        for x, y in product(*(partitions(n) for n in scalars)):
            cycles = tuple(sorted(x + y + (1,)))
            result += Fraction(decomposition('F4', cycles).get(ZERO, 0),
                               centralizer(x) * centralizer(y))
        assert result.denominator == 1
        return int(result)

    tail_weights = [source_weight(counts) for counts in ((3, 1), (4, 0))]
    paths = []
    for partition in partitions(5):
        shape = partition[::-1]
        multiplicity = schur_multiplicity(shape)
        if multiplicity:
            for pure in horizontal_strip(shape, 4):
                for marked in horizontal_strip(pure, 3):
                    if transpose((4, 4, 4, 2)) in horizontal_strip(transpose(marked), 2):
                        paths.append([shape, pure, marked, multiplicity])
    assert tail_weights == [2, 0]
    assert paths == [[(3, 1, 1), (4, 3, 1, 1), (4, 4, 3, 1), 2]]
    return dict(tail_source_weights=tail_weights,
                five_letter_shapes={str(s[::-1]): schur_multiplicity(s[::-1])
                                    for s in partitions(5) if schur_multiplicity(s[::-1])},
                arm_paths=paths, source_upper_bounds=[2, 2])


def ambient(archive):
    shape = (4, 4, 4, 2)
    def singlet(cycles):
        return decomposition('F4', cycles).get(ZERO, 0)
    sizes = [multiplicity(transpose(shape), n, singlet) for n in (12, 13, 14)]
    # Check the terminal and one-double-jet terms directly in odd Taylor
    # directions, independently of the symmetric-jet character calculation.
    source = 0
    for parts in partitions(12):
        smaller = parts[::-1]
        if transpose(shape) in horizontal_strip(transpose(smaller), 2):
            source += sum(Fraction(character(smaller, cycles) * singlet((1,) + cycles),
                                   centralizer(cycles)) for cycles in partitions(12))
    assert sizes == [10007, 3720, 704]
    assert source == sizes[1] and schur_multiplicity(shape) == sizes[2]
    matrices, _ = archive_matrices(archive)
    cartan = matrices[:4]
    assert all(np.array_equal(m, np.diag(m.diagonal())) for m in cartan)
    form = np.einsum('aij,bji->ab', cartan, cartan)
    assert not (form % 6).any()
    form //= 6
    four = int(fmpz_mat(form.tolist()).det())
    three = int(fmpz_mat(form[:3,:3].tolist()).det())
    assert (four, three) == (4, 8)
    return dict(lengths=[12, 13, 14], dimensions=sizes, cartan_gram=[four, three],
                cartan_product=four*three, outgoing_rank_upper=sizes[-1]-1)


def sign(sequence):
    return (-1) ** sum(a > b for i, a in enumerate(sequence) for b in sequence[i+1:])


def crossed_terms():
    # Crossed terms in the polarized quartic contain Tr([A,B]C), a
    # three-form. Keep each symmetric Jordan arm and each free second jet
    # formal, so cancellation here is an integer polynomial identity.
    closed, marked = Counter(), Counter()
    for p, q, r, s in product(permutations(range(4)), permutations(range(4)),
                              permutations(range(3)), permutations(range(3))):
        arms = tuple(sorted(tuple(sorted(row)) for row in zip(p[:2],q[:2],r[:2],s[:2])))
        weight = sign(p)*sign(q)*sign(r)*sign(s)
        for c, d in ((r[2], s[2]), (s[2], r[2])):
            v = q[2:] + (d,)
            if len(set(v)) != 3:
                continue
            weighted = weight * sign(v)
            marked[arms, tuple(sorted(p[2:])), c, tuple(sorted(v))] += weighted * sign(p[2:])
            u = p[2:] + (c,)
            if len(set(u)) == 3:
                forms = tuple(sorted((tuple(sorted(u)), tuple(sorted(v)))))
                closed[arms, forms] += weighted * sign(u)
    assert not any(closed.values()) and not any(marked.values())
    return dict(tree_monomials_tested=len(closed), source_monomials_tested=len(marked),
                nonzero_residuals=[0, 0])


TRIPLES = tuple(combinations_with_replacement(range(4), 3))
QUARTICS = tuple(combinations_with_replacement(range(4), 4))
BRACKETS = tuple(combinations(range(4), 2))
TAILS = tuple(product(range(4), repeat=5))


@cache
def marked_layout():
    # Each tuple is one column of the diagram. The pure arm is the horizontal
    # four-strip, the three ordinary slots on the other arm the three-strip,
    # and the two slots of the doubled jet the final vertical strip.
    groups = Counter()
    for p, q, r, s in product(permutations(range(4)), permutations(range(4)),
                              permutations(range(3)), permutations(range(3))):
        # The vertical strip occupies (row 2,column 3) and (row 3,column 1)
        # in zero-based order.
        pair = s[2], q[3]
        if pair[0] == pair[1]:
            continue
        pure = tuple(sorted((s[0], q[1], r[1], p[3])))
        ordinary = tuple(sorted((s[1], q[2], r[2])))
        tail = p[0], q[0], r[0], p[1], p[2]
        key = (QUARTICS.index(pure), TRIPLES.index(ordinary), BRACKETS.index(tuple(sorted(pair))),
               TAILS.index(tail))
        groups[key] += sign(p)*sign(q)*sign(r)*sign(s)*sign(pair)
    return np.asarray([(*k, v) for k, v in groups.items() if v], dtype=np.int64).T


@cache
def tail_trace_indices(words):
    return tuple(tuple(np.array([sum(t[i] * 4**power for power, i in enumerate(w[::-1]))
                                 for t in TAILS]) for w in traces) for traces in words)


def marked_coefficients(colours, jordan, inverse, prime, words):
    pairs = (colours[:, None] @ colours[None, :]) % prime
    dim = colours.shape[-1]
    pairs = pairs.reshape(16, dim, dim)
    triples = ((pairs[:, None] @ colours[None, :]) % prime).reshape(64, dim, dim)
    shorter = [pairs, triples]
    trace_values = {k: np.trace(a, axis1=1, axis2=2) % prime for k, a in zip((2, 3), shorter)}
    trace_values[5] = (np.einsum('aij,bji->ab', pairs, triples) % prime).ravel()
    pure_words = np.stack([pairs[4*a+b] @ pairs[4*c+d] % prime for a,b,c,d in QUARTICS])
    pure = np.einsum('tij,bji->bt', jordan, pure_words) % prime
    brackets = np.stack([(pairs[4*a+b] - pairs[4*b+a]) % prime for a,b in BRACKETS])
    triple_words = triples[[16*a + 4*b + c for a,b,c in TRIPLES]]
    marked = (triple_words[:, None] @ brackets[None, :]) % prime
    projected = np.einsum('tij,bcji->bct', jordan, marked) % prime
    pairings = ((pure * inverse % prime) @ projected.reshape(-1, len(jordan)).T) % prime
    pairings = pairings.reshape(len(QUARTICS), len(TRIPLES), len(BRACKETS))
    a, b, bracket, tail, signs = marked_layout()
    weights = np.zeros(len(TAILS), dtype=np.int64)
    np.add.at(weights, tail, pairings[a, b, bracket] * signs)
    weights %= prime
    out = []
    for traces, indices in zip(words, tail_trace_indices(words), strict=True):
        values = np.ones(len(TAILS), dtype=np.int64)
        for word, ix in zip(traces, indices, strict=True):
            values = (values * trace_values[len(word)][ix]) % prime
        out.append(int(weights @ values % prime))
    assert max(dim**2, len(TAILS)) * (prime-1)**2 < 2**63
    return out


def literal_marked(colours, jordan, inverse, prime, words):
    def boxes(shape):
        return {(i,j) for i, length in enumerate(shape) for j in range(length)}
    diagrams = [boxes(shape) for shape in ((3,1,1), (4,3,1,1), (4,4,3,1), (4,4,4,2))]
    parts = [sorted(diagrams[0])] + [sorted(b-a) for a,b in zip(diagrams, diagrams[1:])]

    @cache
    def matrix(labels):
        value = colours[labels[0]]
        for label in labels[1:]:
            value = (value @ colours[label]) % prime
        return value

    @cache
    def pairing(pure, ordinary, double):
        commutator = matrix(double) - matrix(double[::-1])
        marked = (matrix(ordinary) @ commutator) % prime
        a = np.einsum('tij,ji->t', jordan, matrix(pure)) % prime
        b = np.einsum('tij,ji->t', jordan, marked) % prime
        return int(((a * inverse % prime) * b).sum() % prime)

    @cache
    def spectator(labels):
        values = []
        for traces in words:
            value = 1
            for word in traces:
                letters = tuple(labels[i] for i in word)
                value = value * int(np.trace(matrix(letters))) % prime
            values.append(value)
        return values

    result = [0] * len(words)
    choices = [tuple(permutations(range(n))) for n in (4,4,3,3)]
    for tableau in product(*choices):
        tail, pure, ordinary, double = [tuple(tableau[j][i] for i,j in cells) for cells in parts]
        if double[0] == double[1]:
            continue
        factor = pairing(pure, ordinary, double)
        for column in tableau:
            factor *= sign(column)
        for i, value in enumerate(spectator(tail)):
            result[i] += factor * value
    return [value % prime for value in result]


def run(archive, prime, count=72):
    assert sha256(archive.read_bytes()).hexdigest() == ARCHIVE
    character_result = dimensions()
    representation, _ = basis()
    jordan = jordan_colours()
    gram = np.einsum('iab,jba->ij', jordan, jordan)
    assert np.array_equal(gram, np.diag(gram.diagonal()))
    inverse = np.array([pow(int(v), -1, prime) for v in gram.diagonal()], dtype=np.int64)
    # Every quartic arm is symmetric. The quadratic projections vanish, and
    # the symmetrized cubic has zero trace by metric transposition.
    assert not np.einsum('iab,jba->ij', jordan, representation).any()
    for matrix in representation:
        assert not np.einsum('iab,jba->ij', jordan, matrix @ representation).any()
    words = tuple(trace_words(5))
    # The native finite colour dimensions in degrees four and five agree
    # with trace evaluation and independent polarized quartic relations.
    free = {row['free'] for row in json.loads((HERE/'dressed-tree.json').read_bytes())['kernel']}
    with ZipFile(archive) as source:
        old = [poly for i, poly in enumerate(columns(source)) if i not in free]
        original = evaluation_rows(source, old, prime, count, representation)
    rows = []
    for previous, sample in zip(original, colour_samples(representation, prime, count), strict=True):
        outside = coefficients(sample, jordan, inverse, prime)
        assert (outside[0] - 24*outside[1] - 12*outside[2]) % prime == 0
        inside = marked_coefficients(sample, jordan, inverse, prime, words)
        if not rows:
            assert inside == literal_marked(sample, jordan, inverse, prime, words)
        rows.append(previous + outside[:2] + inside)
    matrix = nmod_mat(rows, prime)
    column_ids = pivots(matrix)
    row_ids = pivots(matrix.transpose())
    assert column_ids[:63] == list(range(63))
    minor = nmod_mat([[rows[i][j] for j in column_ids] for i in row_ids], prime)
    determinant = int(minor.det())
    assert determinant and len(column_ids) == 61 + sum(character_result['source_upper_bounds'])
    return dict(prime=prime, archive_sha256=ARCHIVE, dimensions=character_result, samples=count,
                old_rank=61, tail_adds=2, joint_rank=len(column_ids),
                ambient=ambient(archive), quartic_crossed_terms=crossed_terms(),
                literal_marked_traces=len(words),
                marked_pivots=[words[i-63] for i in column_ids[63:]], minor_columns=column_ids,
                minor_rows=row_ids, determinant=determinant)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--prime', type=int, required=True, choices=(1000003, 1000033))
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    print(json.dumps(run(args.archive, args.prime), separators=(',', ':')))
