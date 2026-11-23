"""Use the Albert multiplication operators in the shortened eleven-letter target."""

import argparse
from functools import cache
import json
from pathlib import Path

from flint import nmod_mat
import numpy as np
from scipy.sparse import csr_matrix

from glueball_factor_test import COLOURS, HEIGHTS, TraceTests, flags, layout
from jordan_bracket_test import albert_basis, jordan_colours
from source_symmetry import multiplicity, singlet_trace


def arm_words(slots, digits):
    word = np.zeros(digits.shape[1], dtype=np.int64)
    reverse = word.copy()
    for i, j in zip(slots, reversed(slots)):
        word = COLOURS * word + digits[i]
        reverse = COLOURS * reverse + digits[j]
    phase = np.sign(reverse - word) if len(slots) % 2 else 1
    return np.minimum(word, reverse), phase


class JordanTests:
    def __init__(self, seeds):
        digits, signs = flags()
        rows = []
        for seed in seeds:
            a, b = layout((5, 6), seed)
            x, sx = arm_words(a, digits)
            y, sy = arm_words(b, digits)
            keys, coordinates = np.unique(COLOURS**6 * x + y, return_inverse=True)
            coeff = np.bincount(coordinates, weights=signs * sx * sy).astype(np.int64)
            rows.append((keys[coeff != 0], coeff[coeff != 0]))
        keys = np.unique(np.concatenate([k for k, _ in rows]))
        offsets = np.cumsum([0] + [len(k) for k, _ in rows])
        cols = np.concatenate([np.searchsorted(keys, k) for k, _ in rows])
        self.coefficients = csr_matrix(
            (np.concatenate([v for _, v in rows]), cols, offsets),
            shape=(len(seeds), len(keys)),
        )
        x, self.left = np.unique(keys // COLOURS**6, return_inverse=True)
        y, self.right = np.unique(keys % COLOURS**6, return_inverse=True)
        self.five = x
        self.six = y
        self.operators = jordan_colours()
        gram = np.einsum("aij,bji->ab", self.operators, self.operators)
        assert np.array_equal(gram, np.diag(gram.diagonal()))
        self.gram = gram.diagonal()

    def point(self, matrices, prime):
        size = matrices.shape[1]
        form = np.array([pow(int(x), -1, prime) for x in self.gram])
        # Tr(L W) uses only the nonzero entries of the multiplication operators.
        projection = csr_matrix(self.operators.swapaxes(1, 2).reshape(len(form), -1))
        five = np.union1d(self.five, self.six // COLOURS)
        words = {5: five}
        for n in range(4, 0, -1):
            words[n] = np.unique(words[n + 1] // COLOURS)
        matrices = matrices.astype(np.float64)
        previous = matrices[words[1]]

        def multiply(left, right):
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                out = left.reshape(-1, size) @ right
            assert 0 <= out.min() and out.max() <= size * (prime - 1) ** 2 < 2**53
            return (out % prime).reshape(-1, size, size)

        for n in range(2, 6):
            w = words[n]
            current = np.empty((len(w), size, size))
            for j in range(COLOURS):
                indices = np.flatnonzero(w % COLOURS == j)
                if len(indices):
                    current[indices] = multiply(
                        previous[np.searchsorted(words[n - 1], w[indices] // COLOURS)],
                        matrices[j],
                    )
            previous = current

        def project(block):
            return (
                projection @ block.reshape(-1, size * size).astype(np.int64).T
            ).T % prime

        short = project(previous[np.searchsorted(five, self.five)]) * form % prime
        long = np.empty((len(self.six), len(form)), dtype=np.int64)
        for j in range(COLOURS):
            positions = np.flatnonzero(self.six % COLOURS == j)
            for start in range(0, len(positions), 128):
                ids = positions[start : start + 128]
                v = multiply(
                    previous[np.searchsorted(five, self.six[ids] // COLOURS)],
                    matrices[j],
                )
                long[ids] = project(v)
        joined = np.empty(len(self.left), dtype=np.int64)
        for i in range(0, len(joined), 1024):
            part = slice(i, i + 1024)
            joined[part] = (short[self.left[part]] * long[self.right[part]]).sum(
                axis=1
            ) % prime
        return self.coefficients @ joined % prime

    def direct(self, seed, matrices, prime):
        digits, signs = flags()
        a, b = layout((5, 6), seed)

        @cache
        def word(w):
            if len(w) == 1:
                return matrices[w[0]]
            return word(w[:-1]) @ matrices[w[-1]] % prime

        def project(slots):
            words, inverse = np.unique(digits[slots].T, axis=0, return_inverse=True)
            values = [
                np.einsum("aij,ji->a", self.operators, word(tuple(w))) for w in words
            ]
            return np.array(values)[inverse] % prime

        form = np.array([pow(int(x), -1, prime) for x in self.gram])
        return int(signs @ (project(a) * form * project(b)).sum(axis=1)) % prime


def main(args):
    ref = json.loads(args.traces.read_bytes())
    old = [(tuple(s), j) for s, j in ref["recipes"]]
    trace_tests = TraceTests(old)
    seeds = list(range(args.seeds))
    jordan = JordanTests(seeds)
    generators, _ = albert_basis()
    samples = list(range(args.sample_start, args.sample_start + args.points))
    frame = []
    for i, seed in enumerate(samples):
        coefficients = np.random.default_rng(seed).integers(-3, 4, size=(COLOURS, 52))
        matrices = (coefficients @ generators.reshape(52, -1)).reshape(
            COLOURS, 27, 27
        ) % args.prime
        extra = jordan.point(matrices, args.prime)
        if i == 0:
            active = np.flatnonzero(extra)
            if len(active):
                j = int(active[0])
                assert extra[j] == jordan.direct(seeds[j], matrices, args.prime)
        frame.append(trace_tests.point(matrices, args.prime).tolist() + extra.tolist())
    matrix = nmod_mat(frame, args.prime)

    def pivots(a):
        r, rank = a.rref()
        return [next(j for j in range(a.ncols()) if r[i, j]) for i in range(rank)]

    ii, jj = pivots(matrix.transpose()), pivots(matrix)
    dimension = multiplicity(HEIGHTS, 11, singlet_trace)
    assert dimension == 55 and len(ii) <= dimension
    det = int(nmod_mat([[frame[i][j] for j in jj] for i in ii], args.prime).det())
    assert det
    result = dict(
        taylor_shape=HEIGHTS,
        terminal_dimension=dimension,
        prime=args.prime,
        rank_on_defining_traces=nmod_mat(
            [r[: len(old)] for r in frame], args.prime
        ).rank(),
        augmented_rank=len(ii),
        determinant=det,
        jordan_seeds_tested=seeds,
        trace_recipes=[old[j] for j in jj if j < len(old)],
        jordan_pivots=[seeds[j - len(old)] for j in jj if j >= len(old)],
        sample_seeds=[samples[i] for i in ii],
    )
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(result["rank_on_defining_traces"], len(ii), det, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traces", type=Path, required=True)
    p.add_argument("--seeds", type=int, default=32)
    p.add_argument("--points", type=int, default=59)
    p.add_argument("--sample-start", type=int, default=954327)
    p.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if not __debug__ or min(a.seeds, a.points) < 1:
        p.error("enable assertions and use positive sample sizes")
    main(a)
