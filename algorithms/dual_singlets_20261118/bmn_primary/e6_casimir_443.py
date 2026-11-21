"""Two classical boundaries from the split Killing Casimir in the 27 of E6."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations
import json
from pathlib import Path

import numpy as np
from flint import nmod_mat
from scipy.sparse import csr_matrix

from e6_scalar_boundaries import frame, generators, m, trace_boundaries


def alternating_words(slots):
    result = Counter()
    for perm in permutations(range(len(slots))):
        words = [[], [], []]
        for slot, letter in zip(slots, perm):
            words[slot].append(letter)
        sign = (-1)**sum(i>k for j,i in enumerate(perm) for k in perm[j+1:])
        result[tuple(map(tuple, words))] += sign
    return result


def expansion(factors):
    result = {((), (), ()): 1}
    for slots in factors:
        terms = Counter()
        for words, coefficient in result.items():
            for step, sign in alternating_words(slots).items():
                terms[tuple(u+v for u,v in zip(words,step))] += coefficient*sign
        result = {w:c for w,c in terms.items() if c}
    return result


def difference(pair):
    triple = (0, 1, 2)
    result = Counter(expansion((triple, pair, triple, triple)))
    result.subtract(expansion((pair, triple, triple, triple)))
    return {w:c for w,c in result.items() if c}


def projections(prime):
    basis = generators()
    # The sum of root outer products is four times the minuscule weight sum.
    gram = 4*np.einsum("iab,jba->ij", basis, basis)
    inverse = np.asarray(nmod_mat(gram.tolist(), prime).inv().tolist(), dtype=np.int64)
    trace = csr_matrix(basis.transpose(0,2,1).reshape(78,-1))
    dual = csr_matrix((inverse@basis.reshape(78,-1)%prime).T)
    return trace, dual


def split_casimir_values(scalars, polynomials, prime, projections_data, outer_parity):
    """Evaluate twice the specified outer projection of an eleven-scalar trace."""
    if outer_parity not in (1, -1):
        raise ValueError("The diagram involution has eigenvalues 1 and -1")
    trace, dual = projections_data

    @lru_cache(None)
    def word(letters):
        if not letters:
            return np.eye(27, dtype=np.int64)
        return word(letters[:-1])@scalars[letters[-1]]%prime

    @lru_cache(None)
    def projected(letters):
        covector = trace@word(letters).ravel()%prime
        return (dual@covector%prime).reshape(27,27)

    @lru_cache(None)
    def contracted(first, middle, last):
        joined = projected(first)@projected(last)%prime
        return int(np.sum(joined.T*word(middle), dtype=np.int64))%prime

    values = []
    for polynomial in polynomials:
        value = 0
        for (first, middle, last), coefficient in polynomial.items():
            forward = contracted(first, middle, last)
            # Charge conjugation gives eleven minus signs and reversed traces.
            reverse = contracted(last[::-1], middle[::-1], first[::-1])
            value += coefficient*(forward-outer_parity*reverse)
        values.append(value%prime)
    return values


def compute(work, prime, seed):
    work.mkdir(parents=True, exist_ok=True)
    letters = [(0,0,(i,)) for i in range(3)]
    points = 16
    evals = frame(letters, points, prime, seed)
    traces = trace_boundaries((4,4,3), 1, evals, m.library(work,27,prime), prime)
    pairs = ((0,2), (0,1))
    polynomials = [difference(pair) for pair in pairs]
    forms = projections(prime)
    extra = [split_casimir_values([evals.L[letter][1][i,0] for letter in letters],
                                  polynomials, prime, forms, 1) for i in range(points)]
    combined = np.vstack((traces, np.array(extra, dtype=np.int64).T))
    record = dict(prime=prime, seed=seed, points=points, trace_parents=len(traces),
                  weight=[4,4,3], length=11, outer=1, tensor_slot_pairs=pairs)
    for name, matrix in (("trace", traces), ("with_casimir", combined)):
        minor = m.minor(matrix, prime)
        record[name] = dict(minor=minor,
                            values=matrix[np.ix_(minor["rows"],minor["columns"])].tolist())
    assert record["trace"]["minor"]["rank"] == 12
    assert record["with_casimir"]["minor"]["rank"] == 14
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--prime", type=int, choices=(1000003,1000033), required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = compute(args.work, args.prime, args.seed)
    args.out.write_text(json.dumps(data,separators=(",",":"))+"\n")
    print(data["with_casimir"]["minor"]["determinant"])
