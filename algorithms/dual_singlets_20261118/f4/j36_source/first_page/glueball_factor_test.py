"""Test the first image of a shortened F4 source on the (5,4,2) Taylor component."""

import argparse
from functools import cache
from itertools import permutations, product
import json
from pathlib import Path

from flint import nmod_mat
import numpy as np
from scipy.sparse import csr_matrix

from flag_boundaries import adjoint
from single_trace_flag import native26
from source_symmetry import multiplicity, partitions, singlet_trace


HEIGHTS = (5, 4, 2)
COLOURS = max(HEIGHTS)


@cache
def flags():
    groups = []
    for height in HEIGHTS:
        signed = []
        for p in permutations(range(height)):
            inversions = sum(x > y for i, x in enumerate(p) for y in p[i + 1 :])
            signed.append((p, (-1) ** inversions))
        groups.append(signed)
    digits, signs = [], []
    for a, b, c in product(*groups):
        digits.append(a[0] + b[0] + c[0])
        signs.append(a[1] * b[1] * c[1])
    assert len(signs) == 5760
    return np.asarray(digits, dtype=np.int64).T.copy(), np.asarray(
        signs, dtype=np.int64
    )


def layout(shape, seed):
    slots = np.random.default_rng(382851 + 1171 * seed).permutation(sum(shape))
    return np.split(slots, np.cumsum(shape)[:-1])


def canonical(slots, digits):
    length = len(slots)
    forward = np.zeros(digits.shape[1], dtype=np.int64)
    reverse = forward.copy()
    for i, j in zip(slots, reversed(slots)):
        forward = COLOURS * forward + digits[i]
        reverse = COLOURS * reverse + digits[j]
    past = COLOURS ** (length - 1)

    def smallest(word):
        best = word.copy()
        for _ in range(length - 1):
            word = (word % past) * COLOURS + word // past
            best = np.minimum(best, word)
        return best

    a, b = smallest(forward), smallest(reverse)
    return np.minimum(a, b), np.sign(b - a) if length % 2 else 1


def terms(shape, seed):
    digits, signs = flags()
    weight = signs.copy()
    codes = []
    for group in layout(shape, seed):
        code, phase = canonical(group, digits)
        weight *= phase
        codes.append(code)
    packed = np.zeros(len(weight), dtype=np.int64)
    for length in dict.fromkeys(shape):
        block = np.sort([c for c, n in zip(codes, shape) if n == length], axis=0)
        for code in block:
            packed = packed * COLOURS**length + code
    keys, columns = np.unique(packed, return_inverse=True)
    coefficients = np.bincount(columns, weights=weight).astype(np.int64)
    keep = coefficients != 0
    return keys[keep], coefficients[keep]


class TraceTests:
    def __init__(self, recipes):
        self.blocks = []
        required = {}
        for shape in sorted({s for s, _ in recipes}):
            indices = [i for i, (s, _) in enumerate(recipes) if s == shape]
            polynomials = [terms(shape, recipes[i][1]) for i in indices]
            union = np.unique(np.concatenate([k for k, _ in polynomials]))
            if not len(union):
                continue
            cols = np.concatenate([np.searchsorted(union, k) for k, _ in polynomials])
            coefficients = np.concatenate([v for _, v in polynomials])
            starts = np.cumsum([0] + [len(k) for k, _ in polynomials])
            matrix = csr_matrix(
                (coefficients, cols, starts), shape=(len(indices), len(union))
            )
            words = []
            for n in reversed(shape):
                factor, inverse = np.unique(union % COLOURS**n, return_inverse=True)
                words.append((n, factor, inverse))
                required.setdefault(n, []).append(factor)
                union = union // COLOURS**n
            self.blocks.append((indices, matrix, words))
        self.words = {
            n: np.unique(np.concatenate(arrays)) for n, arrays in required.items()
        }
        self.size = len(recipes)
        for _, _, words in self.blocks:
            for i, (n, w, where) in enumerate(words):
                words[i] = (n, np.searchsorted(self.words[n], w)[where])

        # Share both halves of a trace across the entire flag basis.
        needed = {}
        for n, words in self.words.items():
            right = n - n // 2
            needed.setdefault(n // 2, []).append(words // COLOURS**right)
            needed.setdefault(right, []).append(words % COLOURS**right)
        self.halves = {}
        self.prefix = {}
        for n in range(max(needed), 0, -1):
            w = np.unique(np.concatenate(needed[n]))
            self.halves[n] = w
            if n > 1:
                needed.setdefault(n - 1, []).append(w // COLOURS)
        for n, w in self.halves.items():
            if n > 1:
                self.prefix[n] = np.searchsorted(self.halves[n - 1], w // COLOURS)
        self.trace_pairs = {}
        for n, w in self.words.items():
            left, right = n // 2, n - n // 2
            self.trace_pairs[n] = (
                np.searchsorted(self.halves[left], w // COLOURS**right),
                np.searchsorted(self.halves[right], w % COLOURS**right),
            )

    def point(self, matrices, prime):
        size = matrices.shape[1]
        assert size**2 * (prime - 1) ** 2 < 2**53
        values = {1: matrices[self.halves[1]].astype(np.float64)}
        for n in range(2, max(self.halves) + 1):
            words = self.halves[n]
            table = np.empty((len(words), size, size), dtype=np.float64)
            for i in range(COLOURS):
                positions = np.flatnonzero(words % COLOURS == i)
                if not positions.size:
                    continue
                with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                    vals = (
                        values[n - 1][self.prefix[n][positions]].reshape(-1, size)
                        @ matrices[i]
                    )
                assert 0 <= vals.min() and vals.max() <= size * (prime - 1) ** 2
                table[positions] = (vals % prime).reshape(-1, size, size)
            values[n] = table
        traces = {}
        for n, (a, b) in self.trace_pairs.items():
            left, right = n // 2, n - n // 2
            entries = np.empty(len(a), dtype=np.int64)
            for i in range(0, len(a), 128):
                block = slice(i, i + 128)
                entries[block] = np.einsum(
                    "ijk,ikj->i", values[left][a[block]], values[right][b[block]]
                )
            traces[n] = entries % prime
        result = np.zeros(self.size, dtype=np.int64)
        for rows, coefficients, factors in self.blocks:
            value = np.ones(coefficients.shape[1], dtype=np.int64)
            for n, columns in factors:
                value = value * traces[n][columns] % prime
            result[rows] = (coefficients @ value) % prime
        return result


def direct(shape, seed, matrices, prime):
    digits, signs = flags()

    @cache
    def matrix(word):
        if len(word) == 1:
            return matrices[word[0]]
        return matrix(word[:-1]) @ matrices[word[-1]] % prime

    weights = signs.copy()
    for group in layout(shape, seed):
        words, inverse = np.unique(digits[group].T, axis=0, return_inverse=True)
        cut = len(group) // 2
        values = [
            int((matrix(tuple(w[:cut])) * matrix(tuple(w[cut:])).T).sum()) % prime
            for w in words
        ]
        weights = weights * np.asarray(values)[inverse] % prime
    return int(weights.sum()) % prime


def main(args):
    if args.replay:
        previous = json.loads(args.replay.read_text())
        recipes = [(tuple(s), j) for s, j in previous["recipes"]]
        samples = previous["sample_seeds"]
    else:
        shapes = [tuple(reversed(p)) for p in partitions(11) if min(p) >= 2]
        recipes = [(s, j) for s in shapes for j in range(args.seeds)]
        samples = list(range(args.sample_start, args.sample_start + args.points))
    dimension = multiplicity(HEIGHTS, 11, singlet_trace)
    assert dimension == 55
    rep = native26() if args.rep == 26 else adjoint()
    tests = TraceTests(recipes)
    frame = []
    for i, seed in enumerate(samples):
        coeffs = np.random.default_rng(seed).integers(-3, 4, size=(COLOURS, 52))
        matrices = (coeffs @ rep.reshape(52, -1)).reshape(
            COLOURS, args.rep, args.rep
        ) % args.prime
        row = tests.point(matrices, args.prime)
        if i == 0:
            active = np.flatnonzero(row)
            for j in active[:: max(1, len(active) // 3)][:3]:
                assert row[j] == direct(*recipes[j], matrices, args.prime)
        frame.append(row.tolist())
    evaluated = nmod_mat(frame, args.prime)

    def pivots(m):
        reduced, rank = m.rref()
        return [next(j for j in range(m.ncols()) if reduced[i, j]) for i in range(rank)]

    rows, cols = pivots(evaluated.transpose()), pivots(evaluated)
    det = int(nmod_mat([[frame[i][j] for j in cols] for i in rows], args.prime).det())
    assert det and len(cols) <= dimension
    result = dict(
        taylor_shape=HEIGHTS,
        scalar_shape=(3, 3, 2, 2, 1),
        terminal_dimension=dimension,
        representation=args.rep,
        prime=args.prime,
        rank=len(cols),
        determinant=det,
        recipes=[recipes[j] for j in cols],
        sample_seeds=[samples[i] for i in rows],
    )
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{len(cols)} of {dimension}, determinant {det}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--rep", type=int, choices=(26, 52), default=26)
    parser.add_argument("--seeds", type=int, default=32)
    parser.add_argument("--points", type=int, default=60)
    parser.add_argument("--sample-start", type=int, default=673201)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or min(args.points, args.seeds) < 1:
        parser.error("enable assertions and use positive dimensions")
    main(args)
