"""Verify the rank of the projected F4 tree images for four scalar dressings."""
import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import permutations, product
import json
from math import prod
from pathlib import Path
from random import Random
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flint import nmod_mat
import numpy as np
from scipy.sparse import csr_matrix

from verify import ARCHIVE, HERE, contents


def nonzero(row):
    return {key: value for key, value in row.items() if value}


@cache
def colour_trace(word):
    forward = min(word[i:] + word[:i] for i in range(len(word)))
    rev = word[::-1]
    backward = min(rev[i:] + rev[:i] for i in range(len(rev)))
    if forward == backward and len(word) % 2:
        return forward, 0
    return (forward, 1) if forward <= backward else (backward, (-1)**len(word))


def put(row, traces, value):
    words, signs = zip(*(colour_trace(trace) for trace in traces))
    row[tuple(sorted(words))] += value * prod(signs)


def multiply(left, right):
    output = Counter()
    for words, a in left.items():
        for more, b in right.items():
            output[tuple(sorted(words + more))] += a*b
    return nonzero(output)


def add_scaled(row, other, factor):
    for k, v in other.items():
        row[k] += factor*v


def bracket(left, right):
    result = Counter()
    for a, x in left.items():
        for b, y in right.items():
            result[a+b] += x*y
            result[b+a] -= x*y
    return nonzero(result)


def trace(*arguments):
    result = Counter()
    for pieces in product(*(arg.items() for arg in arguments)):
        words, coefficients = zip(*pieces)
        put(result, (sum(words, ()),), prod(coefficients))
    return nonzero(result)


def dressings():
    scalars = [{(i,): 1} for i in range(3)]
    commutators = [bracket(scalars[(i+1)%3], scalars[(i+2)%3]) for i in range(3)]
    gram = Counter()
    for perm in permutations(range(3)):
        sign = (-1)**sum(a>b for i,a in enumerate(perm) for b in perm[i+1:])
        term = {(): 1}
        for i,j in enumerate(perm):
            term = multiply(term, trace(scalars[i],scalars[j]))
        add_scaled(gram,term,sign)
    cubic = trace(scalars[0],commutators[0])
    mixed = Counter()
    for i,j in product(range(3),repeat=2):
        add_scaled(mixed, multiply(trace(scalars[i],scalars[j]),
                                   trace(commutators[i],commutators[j])), 1)
    return (nonzero(gram), multiply(cubic,cubic), nonzero(mixed),
            trace(commutators[0],bracket(commutators[1],commutators[2])))


@cache
def interchanges(words):
    marked, scalar = [], []
    for i, word in enumerate(words):
        for j, letter in enumerate(word):
            (marked if letter == 3 else scalar).append((i,j))
    assert len(marked) == 2 and len(scalar) == 12
    result = Counter()
    for (i,j),(k,l) in product(marked,scalar):
        swapped = [list(t) for t in words]
        swapped[i][j],swapped[k][l] = swapped[k][l],swapped[i][j]
        put(result, tuple(map(tuple,swapped)), 1)
    return nonzero(result)


def projected(row):
    # The ordinary GL(4) transfer has diagonal term 6.  Multiply the
    # projector by 84 so both the source coefficients and all minors are integral.
    first, second = Counter(), Counter()
    for w,c in row.items():
        add_scaled(first,interchanges(w),c)
    for w,c in nonzero(first).items():
        add_scaled(second,interchanges(w),c)
    add_scaled(second,first,-8)
    return nonzero(second)


def columns(archive):
    bank = json.loads(contents(archive,'cochains/project-bank.json'))['quotient_seeds']
    factors = dressings()
    for src in bank:
        output = Counter()
        for traces,c in src['out']:
            put(output, tuple(tuple(min(k,3) for k in t) for t in traces), int(c))
        for factor in factors:
            yield projected(multiply(factor,nonzero(output)))
        interchanges.cache_clear()


def colour_samples(representation, prime, count):
    rng = Random(9311422901)
    colours, dim, _ = representation.shape
    for _ in range(count):
        coords = np.array([[rng.randrange(-4,5) for _ in range(colours)] for _ in range(4)])
        yield (coords @ representation.reshape(colours,-1)).reshape(4,dim,dim) % prime


def evaluation_rows(archive, polynomials, prime, count, representation=None):
    if representation is None:
        representation = np.zeros((52,26,26),dtype=np.int64)
        matrices = json.loads(contents(archive,'collected-O0/color-basis-native.json'))['representation26']
        for i,entries in enumerate(matrices):
            for a,b,c in entries:
                representation[i,a,b] = int(c)
    monomials = sorted(set().union(*(p.keys() for p in polynomials)))
    index = {word:i for i,word in enumerate(monomials)}
    data, indices, offsets = [], [], [0]
    for c in polynomials:
        for word,n in c.items():
            data.append(n%prime)
            indices.append(index[word])
        offsets.append(len(data))
    coefficients = csr_matrix((data,indices,offsets),shape=(len(polynomials),len(monomials)),
                               dtype=np.int64)
    trace_words = sorted({trace for word in monomials for trace in word})
    positions = {t:i for i,t in enumerate(trace_words)}
    lookup = np.full((max(map(len,monomials)),len(monomials)),len(trace_words),dtype=int)
    for j,word in enumerate(monomials):
        for i,t in enumerate(word):
            lookup[i,j] = positions[t]
    assert int(coefficients.sum(axis=1).max())*(prime-1) < 2**63
    assert representation.shape[1]**2*(prime-1)**2 < 2**63
    output = []
    for colours in colour_samples(representation,prime,count):

        @cache
        def matrix(word):
            if len(word)==1:
                return colours[word[0]]
            return (matrix(word[:-1]) @ colours[word[-1]])%prime

        def evaluate_trace(word):
            h = len(word)//2
            return int(np.sum(matrix(word[:h])*matrix(word[h:]).T))%prime

        values = np.array([evaluate_trace(t) for t in trace_words]+[1],dtype=np.int64)
        matrix.cache_clear()
        row = np.ones(len(monomials),dtype=np.int64)
        for ix in lookup:
            row = row*values[ix]%prime
        output.append(((coefficients @ row)%prime).tolist())
    return output


def verify(archive, certificate):
    assert certificate['archive_sha256'] == ARCHIVE
    polynomials = list(columns(archive))
    assert len(polynomials) == 80
    basis = certificate['kernel']
    free = [s['free'] for s in basis]
    assert len(set(free)) == len(free) == 19
    for vector in basis:
        coefficients = dict(vector['coefficients'])
        assert coefficients[vector['free']] and set(coefficients).intersection(free)=={vector['free']}
    pivots = [i for i in range(80) if i not in free]
    for record in certificate['minors']:
        p = record['prime']
        rows = evaluation_rows(archive, polynomials, p, max(record['rows'])+1)
        values = [[rows[i][j] for j in pivots] for i in record['rows']]
        assert len(values) == len(pivots) == 61
        assert int(nmod_mat(values,p).det()) == record['determinant'] != 0

    # The pinned modules reduce trace words by the polarized quartic relation.
    # check_colour derives that relation from the supplied 26-dimensional matrices.
    from native import check_colour
    assert check_colour(json.loads(contents(archive,'collected-O0/color-basis-native.json')))[0]['cartan_quartic_coefficients']==35
    with TemporaryDirectory() as temp:
        for name in ('trace_normal_form','quartic_library'):
            (Path(temp)/(name+'.py')).write_bytes(contents(archive,'cochains/source-replay/'+name+'.py'))
        sys.path.insert(0,temp)
        try:
            from quartic_library import red
            for vector in basis:
                formal = Counter()
                for k,v in vector['coefficients']:
                    add_scaled(formal,polynomials[k],v)
                formal = nonzero(formal)
                assert len(formal) == vector['formal_terms']
                result = Counter()
                for word,v in formal.items():
                    for pieces in product(*(red(t).items() for t in word)):
                        traces, factors = zip(*pieces)
                        result[tuple(sorted(sum(traces,())))] += v*prod(factors)
                assert not nonzero(result)
        finally:
            sys.path.remove(temp)
    return len(pivots), len(basis)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    args=parser.parse_args()
    if not __debug__ or sha256(args.archive.read_bytes()).hexdigest()!=ARCHIVE:
        parser.error('use the pinned archive and enabled assertions')
    with ZipFile(args.archive) as z:
        result=verify(z,json.loads((HERE/'dressed-tree.json').read_bytes()))
    print('Projected tree rank and kernel dimension over Q:',*result)
