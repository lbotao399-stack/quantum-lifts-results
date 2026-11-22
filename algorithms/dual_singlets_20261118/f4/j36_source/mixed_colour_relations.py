"""Compare polarized quartics with the full F4 invariant dimension at low degree."""
import argparse
from collections import Counter
from functools import cache
from itertools import permutations, product
import json
from math import prod
from pathlib import Path
from random import Random
import sys
from zipfile import ZipFile
from hashlib import sha256

from flint import nmod_mat
import numpy as np

from dressed_tree import colour_trace, nonzero
from verify import ARCHIVE, HERE, contents

sys.path.insert(0, str(HERE.parents[1]))
from b2_schur_complex import partitions
from f4_singlet_characters import decomposition


def trace_words(degree):
    words = set()
    for lengths in partitions(degree, 2):
        for order in permutations(range(degree)):
            offset = 0
            traces = []
            for n in lengths:
                traces.append(colour_trace(order[offset:offset+n])[0])
                offset += n
            words.add(tuple(sorted(traces)))
    return sorted(words)


def signed_transpositions(words):
    lookup = {word:i for i,word in enumerate(words)}
    degree = sum(map(len, words[0]))
    for i in range(degree-1):
        indices, signs = [], []
        labels = [i+1 if j == i else i if j == i+1 else j for j in range(degree)]
        for word in words:
            parts = [colour_trace(tuple(labels[j] for j in t)) for t in word]
            indices.append(lookup[tuple(sorted(t for t,_ in parts))])
            signs.append(prod(s for _,s in parts))
        yield indices, np.array(signs, dtype=np.int64)


def bracket_word(word):
    terms = {(word[0],):1}
    for i in word[1:]:
        terms = {w+(i,):v for w,v in terms.items()} | {(i,)+w:-v for w,v in terms.items()}
    return terms


def generators(degree, words):
    for active in range(4,degree+1):
        for lengths in partitions(active):
            if len(lengths) != 4:
                continue
            blocks = []
            offset = 0
            for n in lengths:
                blocks.append(bracket_word(tuple(range(offset,offset+n))))
                offset += n
            for tail in partitions(degree-active,2):
                rest = []
                offset = active
                for n in tail:
                    rest.append(tuple(range(offset,offset+n)))
                    offset += n
                relation = Counter()
                for pieces in product(*(b.items() for b in blocks)):
                    terms, values = zip(*pieces)
                    factor = prod(values)
                    for a,b,c,d in permutations(terms):
                        for traces,k in [([a+b+c+d],12),([a+b,c+d],-1)]:
                            parts = [colour_trace(w) for w in traces+rest]
                            relation[tuple(sorted(t for t,_ in parts))] += k*factor*prod(s for _,s in parts)
                if nonzero(relation):
                    yield [relation[w] for w in words]


def pivots(matrix):
    reduced, rank = matrix.rref()
    chosen = []
    j = 0
    for i in range(rank):
        while not reduced[i,j]:
            j += 1
        chosen.append(j)
        j += 1
    return chosen


def independent(rows, prime):
    """Keep original rows. Pivot columns of the transpose are independent rows."""
    return rows[pivots(nmod_mat(rows.T.tolist(), prime))]


def orbit_rank(initial, action, prime, upper=None):
    if not initial:
        return 0
    selected = independent(np.asarray(initial, dtype=np.int64) % prime, prime)
    start = 0
    while start < len(selected):
        if len(selected) == upper:
            break
        frontier = selected[start:]
        old = len(selected)
        expanded = [selected] + [frontier[:,ids]*signs % prime for ids,signs in action]
        selected = independent(np.vstack(expanded), prime)
        start = old
    return len(selected)


def sample_colours(matrices, degree, sample, prime):
    rng = Random(65782 + degree*157 + sample*8009)
    weights = np.array([[rng.randrange(-4,5) for _ in matrices] for _ in range(degree)],
                       dtype=np.int64)
    dim = matrices.shape[1]
    colours = (weights @ matrices.reshape(len(matrices),-1)).reshape(degree,dim,dim) % prime
    assert dim*(prime-1)**2 < 2**63
    return colours


def evaluate_words(words, colours, prime):
    @cache
    def matrix(word):
        if len(word) == 1:
            return colours[word[0]]
        return (matrix(word[:-1]) @ colours[word[-1]]) % prime

    @cache
    def trace(word):
        cut = len(word)//2
        # Reduce the diagonal before adding, to keep the bound valid.
        return int(((matrix(word[:cut]) * matrix(word[cut:]).T).sum(axis=0) % prime).sum() % prime)

    return [prod(trace(w) for w in word) % prime for word in words]


def native_row(words, matrices, degree, sample, prime):
    return evaluate_words(words, sample_colours(matrices,degree,sample,prime), prime)


def jordan_completion(words, counts, invariant_dimension, prime):
    from f4_colour import basis, jordan_colours
    derivations, _ = basis()
    jordan = jordan_colours()
    gram = np.einsum('iab,jba->ij',jordan,jordan)
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    inverse = np.array([pow(int(v),-1,prime) for v in gram.diagonal()],dtype=np.int64)
    pairs = []
    for distribution in product(*(range(n+1) for n in counts)):
        if sum(distribution) == 4:
            first = tuple(i for i,n in enumerate(distribution) for _ in range(n))
            second = tuple(i for i,n in enumerate(distribution) for _ in range(counts[i]-n))
            if first <= second:
                pairs.append((first,second))
    pairs.sort()
    rows = []
    for sample in range(invariant_dimension+12):
        colours = sample_colours(derivations,len(counts),sample,prime)
        covariants = {}
        for word in {w for pair in pairs for w in pair}:
            matrix = colours[word[0]]
            for i in word[1:]:
                matrix = (matrix @ colours[i]) % prime
            covariants[word] = np.einsum('iab,ba->i',jordan,matrix) % prime
        extra = [int(((covariants[a]*inverse % prime)*covariants[b]).sum() % prime) for a,b in pairs]
        rows.append(evaluate_words(words,colours,prime)+extra)
    joint = nmod_mat(rows,prime)
    columns = pivots(joint)
    chosen = pivots(joint.transpose())
    minor = nmod_mat([[rows[i][j] for j in columns] for i in chosen],prime)
    determinant = int(minor.det())
    native_rank = sum(i < len(words) for i in columns)
    assert determinant
    assert len(columns) <= invariant_dimension
    return dict(rank=len(columns), defining_rank=native_rank,
                added_pairs=[pairs[j-len(words)] for j in columns if j >= len(words)],
                determinant=determinant, selected_rows=chosen, selected_columns=columns)


def two_matrix_weights(labelled, prime):
    from fractions import Fraction
    from b2_schur_complex import centralizer
    from f4_singlet_characters import ZERO
    result = []
    for a,b in ((5,3),(4,4)):
        labels = (0,)*a+(1,)*b
        words = sorted({t for w in labelled for t,s in [identify(w,labels)] if s})
        dimension = sum(Fraction(decomposition('F4', tuple(sorted(x+y))).get(ZERO,0),
                                 centralizer(x)*centralizer(y)) for x in partitions(a) for y in partitions(b))
        assert dimension.denominator == 1
        d = jordan_completion(words,(a,b),int(dimension),prime)
        # Bound the defining-trace rank from above by actual quartic identities.
        assignments = sorted(set(permutations(labels)))
        rank = specialized_quartic_rank(labelled,words,assignments,prime)
        assert d['defining_rank']+rank == len(words) and d['rank'] == dimension
        result.append(dict(weight=[a,b], stable_traces=len(words), quartic_rank_mod_p=rank,
                           invariant_dimension=int(dimension), **d))
    return result


def specialized_quartic_rank(labelled, words, assignments, prime):
    index = {w:i for i,w in enumerate(words)}
    selected = np.empty((0,len(words)), dtype=np.int64)
    for template in generators(8,labelled):
        terms = [(w,c) for w,c in zip(labelled,template) if c]
        group = set()
        for assign in assignments:
            row = Counter()
            for word, value in terms:
                t, s = identify(word, assign)
                if s:
                    row[index[t]] += s*value
            row = tuple(sorted(nonzero(row).items()))
            if row:
                group.add(row)
        if group:
            extra = np.zeros((len(group),len(words)), dtype=np.int64)
            for i,row in enumerate(sorted(group)):
                for j,v in row:
                    extra[i,j] = v % prime
            selected = independent(np.vstack((selected,extra)), prime)
    return len(selected)


def identify(word, labels):
    parts = [colour_trace(tuple(labels[i] for i in t)) for t in word]
    return tuple(sorted(t for t,_ in parts)), prod(s for _,s in parts)


def paired_degree_eight(prime, matrices, adjoints):
    from math import comb
    from f4_singlet_characters import ZERO
    colour_count = sum(comb(4,k)*decomposition('F4', (2,)*k+(1,)*(8-2*k)).get(ZERO,0)
                       for k in range(5))
    invariant_dimension, remainder = divmod(colour_count,16)
    assert remainder == 0
    labelled = trace_words(8)
    labels = (0,0,1,1,2,2,3,3)
    words = sorted({t for w in labelled for t,s in [identify(w,labels)] if s})
    assignments = sorted(set(permutations(labels)))
    rank = specialized_quartic_rank(labelled,words,assignments,prime)
    one_trace = [w for w in words if len(w) == 1]
    extra_words = Random(65824).sample(one_trace, min(48,len(one_trace)))
    rows = [native_row(words,matrices,4,i,prime) + native_row(extra_words,adjoints,4,i,prime)
            for i in range(invariant_dimension+12)]
    joint = nmod_mat(rows,prime)
    columns = pivots(joint)
    native_rank = sum(j < len(words) for j in columns)
    controls = [j-len(words) for j in columns if j >= len(words)]
    chosen_rows = pivots(joint.transpose())
    minor = nmod_mat([[rows[i][j] for j in columns] for i in chosen_rows],prime)
    determinant = int(minor.det())
    assert determinant and len(columns) <= invariant_dimension
    assert rank+native_rank == len(words)
    jordan = jordan_completion(words,(2,2,2,2),invariant_dimension,prime)
    assert jordan['defining_rank'] == native_rank and jordan['rank'] == invariant_dimension
    return dict(slots=8, multiplicities=[2,2,2,2], stable_traces=len(words),
                quartic_rank_mod_p=rank, native_rank_mod_p=native_rank,
                invariant_dimension=invariant_dimension, augmented_rank=len(columns),
                added_adjoint_words=[extra_words[i][0] for i in controls],
                determinant=determinant, selected_rows=chosen_rows, selected_columns=columns,
                jordan=jordan, two_matrix=two_matrix_weights(labelled,prime))


def table(degree, prime, matrices):
    words = trace_words(degree)
    action = list(signed_transpositions(words))
    rows = list(generators(degree, words))
    rank = orbit_rank(rows, action, prime)
    # All irreducible F4 characters are self-dual.
    left = decomposition('F4', (1,)*(degree//2))
    right = decomposition('F4', (1,)*(degree-degree//2))
    invariant_dimension = sum(c*right.get(w,0) for w,c in left.items())
    evaluations = [native_row(words, matrices, degree, i, prime) for i in range(4)]
    native = orbit_rank(evaluations, action, prime, min(len(words)-rank, invariant_dimension))
    return dict(slots=degree, stable_traces=len(words), quartic_templates=len(rows),
                quartic_rank_mod_p=rank, native_rank_mod_p=native,
                invariant_dimension=invariant_dimension)


def representation(path):
    assert sha256(path.read_bytes()).hexdigest() == ARCHIVE
    with ZipFile(path) as archive:
        data = json.loads(contents(archive, 'collected-O0/color-basis-native.json'))
    result = []
    for key, dim in [('representation26',26), ('adjointMatrices',52)]:
        matrices = np.zeros((52,dim,dim), dtype=np.int64)
        for matrix, entries in zip(matrices, data[key], strict=True):
            for i,j,v in entries:
                matrix[i,j] = int(v)
        result.append(matrices)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--through', type=int, choices=range(4,8), default=7)
    parser.add_argument('--paired-eight', action='store_true')
    parser.add_argument('--prime', type=int, choices=(1000003,1000033), default=1000003)
    args = parser.parse_args()
    if not __debug__:
        parser.error('run with assertions enabled')
    matrices, adjoints = representation(args.archive)
    for degree in range(4,args.through+1):
        print(json.dumps(table(degree,args.prime,matrices)), flush=True)
    if args.paired_eight:
        print(json.dumps(paired_degree_eight(args.prime,matrices,adjoints)), flush=True)
