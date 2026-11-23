"""Test the degree-nine F4 tensor of slot type (3,3,3) on a decorated adjoint edge."""
import argparse
from collections import Counter
from functools import cache
from itertools import permutations, product
import json
from math import prod
from pathlib import Path

from flint import fmpz_mat, nmod_mat
import numpy as np

from jordan_bracket_singlets import sign
from jordan_bracket_test import albert_basis, jordan_colours, native
from source_symmetry import centralizer, character, partitions, singlet_trace
from fractions import Fraction


def alternation():
    rows = Counter()
    for columns in product(permutations(range(3)), repeat=3):
        a, b, c = zip(*columns)
        rows[tuple(sorted(a)), b, tuple(sorted(c))] += np.prod([sign(s) for s in columns])
    return {k: int(v) for k, v in rows.items() if v}


def colours_at(representation, sample):
    rng = np.random.default_rng(933768 + sample)
    coeff = rng.integers(-3, 4, size=(3, 52))
    return (coeff @ representation.reshape(52, -1)).reshape(3, 27, 27)


def integer_traces(words, coefficients, colours):
    matrices = [fmpz_mat(a.tolist()) for a in colours]
    @cache
    def matrix(w):
        if len(w)==1:
            return matrices[w[0][2][0]]
        return matrix(w[:-1]) * matrices[w[-1][2][0]]
    @cache
    def trace(w):
        m=matrix(w)
        return int(sum(m[i,i] for i in range(27)))
    values = [prod(trace(t) for t in word) for word in words]
    return [sum(n*values[i] for i,n in row) for row in coefficients]


def decorated(colours, representation, jordan, inverse, jordan_inverse, prime, rows):
    size = 27
    matrices = representation % prime
    flat = matrices.swapaxes(1, 2).reshape(52, -1)
    transpose = jordan.swapaxes(1, 2).reshape(26, -1)
    actions = []
    for x in colours:
        comm = (x @ matrices - matrices @ x) % prime
        action = inverse @ (flat @ comm.reshape(52, -1).T % prime) % prime
        assert np.array_equal(action.T @ matrices.reshape(52, -1) % prime, comm.reshape(52, -1))
        actions.append(action)
    polynomials = {}
    for letters in sorted({u for a, _, c in rows for u in (a, c)}):
        i, j, k = letters
        word = (colours[i] @ colours[j] % prime) @ colours[k] % prime
        polynomials[letters] = transpose @ (word @ matrices % prime).reshape(52, size**2).T % prime
    total = 0
    for (a, (i, j, k), c), v in rows.items():
        edge = ((actions[i] @ actions[j] % prime) @ actions[k] % prime) @ inverse % prime
        right = polynomials[c] * jordan_inverse[:, None] % prime
        total += v * int(((polynomials[a] @ edge % prime)*right).sum() % prime)
    # The inverse adjoint metric B^{-1} is six times the inverse trace metric.
    return 6 * total % prime


def run(prime):
    shape = (3,3,3)
    upper = sum(Fraction(character(shape, p) * singlet_trace(p), centralizer(p)) for p in partitions(9))
    assert upper == 11
    restrictions = []
    for a,b in (((1,1),(2,2,2)),((2,1),(2,2,1))):
        mult = sum(Fraction(character(a,x)*character(b,y)*character(shape,tuple(sorted((*x,*y,1)))),
                            centralizer(x)*centralizer(y)) for x in partitions(sum(a)) for y in partitions(sum(b)))
        assert mult == 1
        restrictions.append(dict(dotted=a,scalar=b,marked_adjoint_slots=1,multiplicity=int(mult)))
    traces = [r for block in native.words_by_shape(9, (0,0,3,3,3)).values() for r in native.singlets(block)]
    basis, _ = albert_basis()
    gram = np.einsum("aij,bji->ab", basis, basis) % prime
    inverse = np.array(nmod_mat(gram.tolist(), prime).inv().tolist(), dtype=np.int64)
    jordan = jordan_colours()
    jgram = np.einsum("aij,aji->a", jordan, jordan)
    ji = np.array([pow(int(v), -1, prime) for v in jgram])
    patterns = alternation()
    data, exact = [], []
    letters = sorted({l for r in traces for word in r for t in word for l in t})
    words, coefficients = native.flatten(traces)
    for sample in range(int(upper) + 5):
        cols = colours_at(basis, sample)
        frame = native.Frame(27, 1, 1, [0], letters, seed=1, p=prime)
        frame.L = {l:(0, (cols[l[2][0]] % prime).reshape(1,1,27,27)) for l in letters}
        values = frame.eval_words(words, dtype=np.int64).ravel()
        row = [sum(n*int(values[i]) for i,n in r) % prime for r in coefficients]
        integral = integer_traces(words, coefficients, cols)
        assert [x % prime for x in integral] == row
        exact.append(integral)
        value = decorated(cols, basis, jordan, inverse, ji, prime, patterns)
        if sample == 0:
            swap = decorated(cols[[1,0,2]], basis, jordan, inverse, ji, prime, patterns)
            sheared = cols.copy()
            sheared[1] += 2*cols[0]
            assert (swap + value) % prime == 0
            assert decorated(sheared, basis, jordan, inverse, ji, prime, patterns) == value
        data.append(row+[value])
    a = nmod_mat(data, prime)
    ids = native.pivotcols(a)
    rows = native.pivotcols(a.transpose())
    selected = nmod_mat([[data[i][j] for j in ids] for i in rows], prime)
    before = nmod_mat([r[:-1] for r in data], prime).rank()
    trace_rank = fmpz_mat(exact).rank()
    assert len(ids) == upper and before == trace_rank
    return dict(prime=prime, slot_shape=shape, native_character=int(upper), physical_restrictions=restrictions,
                formal_trace_singlets=len(traces),
                defining_trace_evaluation_rank=trace_rank, augmented_rank=len(ids), reduced_alternation_terms=len(patterns),
                column_pivots=ids,row_pivots=rows, determinant=int(selected.det()),
                evaluated_edge=[r[-1] for r in data])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prime",type=int,choices=(1009,1013),default=1009)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    result = run(a.prime)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result["defining_trace_evaluation_rank"],result["augmented_rank"],result["native_character"])
