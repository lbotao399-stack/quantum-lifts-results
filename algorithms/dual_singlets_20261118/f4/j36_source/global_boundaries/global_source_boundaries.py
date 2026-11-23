"""Global classical boundaries in the F4 source of the second differential into S^6."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache
import gzip
from itertools import combinations, permutations, product
import json
from multiprocessing import get_context
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HEIGHTS = (4, 4, 3, 3)
PAIRS = tuple(combinations(range(4), 2))
OFFSETS = (0, 4, 8, 11)
FULL = (1 << 14) - 1
# Counts of the four columns on the left of the colour contraction, arm-bank size,
# and whether to allow a trace or an invariant colour projection inside an arm.
FAMILIES = (
    ((4, 0, 3, 0), 80, False),
    ((4, 0, 2, 1), 44, False),
    ((4, 0, 1, 2), 44, True),
    ((4, 0, 3, 0), 84, True),
    ((4, 4, 0, 0), 65, True),
    ((3, 1, 3, 0), 44, True),
    ((3, 1, 3, 0), 44, False),
    ((4, 1, 2, 0), 46, True),
    ((4, 1, 2, 0), 48, False),
    ((4, 2, 0, 1), 48, True),
)


def inverse(a, p):
    size = len(a)
    aug = np.concatenate((a, np.eye(size, dtype=np.int64)), axis=1) % p
    for j in range(size):
        i = j + int(np.flatnonzero(aug[j:, j])[0])
        aug[[i, j]] = aug[[j, i]]
        row = aug[j] * pow(int(aug[j, j]), -1, p) % p
        aug = (aug - aug[:, j, None] * row) % p
        aug[j] = row
    assert np.array_equal(aug[:, :size], np.eye(size))
    return aug[:, size:]


@cache
def colour(p):
    data = json.loads(gzip.decompress((HERE / "representation26.json.gz").read_bytes()))
    g = np.zeros((52, 26, 26), dtype=np.int64)
    for i, entries in enumerate(data["matrices"]):
        for j, k, value in entries:
            g[i, j, k] = value
    h = np.zeros((26, 26), dtype=np.int64)
    for i, j, value in data["form"]:
        q = Fraction(value)
        h[i, j] = q.numerator * pow(q.denominator, -1, p) % p
    cubic = np.zeros((26, 26, 26), dtype=np.int64)
    data = json.loads((HERE / "global_boundary_cubic.json").read_bytes())
    for i, j, k, value in data["symmetric_entries"]:
        for indices in set(permutations((i, j, k))):
            cubic[indices] = value
    for generator in g:
        action = np.einsum("ri,rjk->ijk", generator, cubic)
        assert not np.any(
            action + action.transpose(1, 2, 0) + action.transpose(2, 0, 1)
        )
        assert not np.any((h @ generator + generator.T @ h) % p)
    jordan = (inverse(h, p) @ cubic) % p
    frames = {}
    for name, mats in [("ad", g), ("jo", jordan)]:
        primal = mats.reshape(len(mats), -1) % p
        dual = mats.transpose(0, 2, 1).reshape(len(mats), -1) % p
        gram = inverse(primal @ dual.T, p)
        frames[name] = (primal.astype(float), dual.T.astype(float), gram.astype(float))
    assert not np.any(
        (g.reshape(52, -1) @ jordan.transpose(0, 2, 1).reshape(26, -1).T) % p
    )
    return g, frames


@cache
def subsets(cuts):
    options = [list(combinations(range(size), k)) for size, k in zip(HEIGHTS, cuts)]
    results = []
    for setting in product(*options):
        mask = sum(
            1 << (OFFSETS[c] + i) for c, used in enumerate(setting) for i in used
        )
        inversions = sum(i - j for used in setting for j, i in enumerate(used))
        results.append((mask, FULL ^ mask, (-1) ** inversions))
    return results


@cache
def extensions(token, mask):
    options = [
        [i for i in range(HEIGHTS[c]) if mask >> (OFFSETS[c] + i) & 1] for c in token
    ]
    result = []
    for labels in product(*options):
        if len(labels) == 2 and labels[0] == labels[1]:
            continue
        previous, phase = mask, 1
        for c, i in reversed(tuple(zip(token, labels))):
            previous ^= 1 << (OFFSETS[c] + i)
            same = (previous >> OFFSETS[c]) & ((1 << HEIGHTS[c]) - 1)
            phase *= (-1) ** ((same >> (i + 1)).bit_count())
        if len(token) == 1:
            index = labels[0]
        else:
            phase *= 1 if labels[0] < labels[1] else -1
            index = 4 + PAIRS.index(tuple(sorted(labels)))
        result.append((previous, phase, index))
    return result


@cache
def arm(spec):
    pattern, k, decoration, stop = spec
    plain = tuple((c,) for c in pattern)
    reverse = plain[:k] + (plain[k + 1], plain[k]) + plain[k + 2 :]
    marked = plain[:k] + (pattern[k : k + 2],) + plain[k + 2 :]
    if decoration:
        cut = stop + (k < stop)
        marked = marked[:stop] + (decoration,) + marked[stop:]
        plain = plain[:cut] + (decoration,) + plain[cut:]
        reverse = reverse[:cut] + (decoration,) + reverse[cut:]
    return marked, plain, reverse, -1 if pattern[k] == pattern[k + 1] else 1


def arm_bank(cuts, count, rng, use_decorations):
    sequence = tuple(i for i, k in enumerate(cuts) for _ in range(k))
    result, seen = [], set()
    while len(result) < count:
        word = tuple(int(x) for x in rng.permutation(sequence))
        k = int(rng.integers(len(word) - 1))
        if word[k] > word[k + 1]:
            word = word[:k] + (word[k + 1], word[k]) + word[k + 2 :]
        decoration, stop = "", 0
        if use_decorations:
            n = int(rng.integers(12))
            if n >= 4:
                decoration = "tr" if n < 7 else ("ad" if n < 9 else "jo")
                stop = int(rng.integers(2, len(word) - 2))
        item = word, k, decoration, stop
        if item not in seen:
            seen.add(item)
            result.append(item)
    return tuple(result)


@cache
def layout(case):
    rng = np.random.default_rng(74935 + case)
    cuts, count, decorations = FAMILIES[case]
    other = tuple(h - n for h, n in zip(HEIGHTS, cuts))
    return (
        cuts,
        arm_bank(cuts, count, rng, decorations),
        arm_bank(other, count, rng, decorations),
    )


class Evaluation:
    def __init__(self, seed, prime, brackets=False):
        generators, self.frames = colour(prime)
        self.p = prime
        coefficients = np.random.default_rng(seed).integers(-3, 4, size=(10, 52))
        self.g = (
            (coefficients @ generators.reshape(52, -1)).reshape(10, 26, 26) % prime
        ).astype(float)
        if brackets:
            for k, (i, j) in enumerate(PAIRS, 4):
                self.g[k] = (self.g[i] @ self.g[j] - self.g[j] @ self.g[i]) % prime
        self.eye = np.eye(26)
        self.zero = self.eye * 0
        # Every GEMM is a sum of integers bounded by 676*(p-1)^2,
        # or by 26*16*(p-1)^2 for the alternations. Reduce between products.
        assert 676 * 24 * prime**2 < 2**53

    def coordinates(self, matrix, name):
        return matrix @ self.frames[name][1] % self.p

    def decorate(self, x, name):
        if name == "tr":
            return self.eye * (x.trace() % self.p)
        primal, _, gram = self.frames[name]
        return (
            ((self.coordinates(x.reshape(676), name) @ gram % self.p) @ primal) % self.p
        ).reshape(26, 26)

    def evaluator(self):
        @cache
        def value(word, mask):
            if not word:
                assert mask == 0
                return self.eye
            token = word[-1]
            if isinstance(token, str):
                return self.decorate(value(word[:-1], mask), token)
            terms = extensions(token, mask)
            if not terms:
                return self.zero
            left = np.concatenate([value(word[:-1], a) for a, _, _ in terms], axis=1)
            right = np.concatenate([self.g[i] * s for _, s, i in terms])
            return (left @ right) % self.p

        return value

    def block(self, case, check=False):
        cuts, left, right = layout(case)
        value = self.evaluator()

        def words(bank, mask):
            f, d = [], []
            for spec in bank:
                marked, plain, reverse, sg = arm(spec)
                f.append(value(marked, mask))
                d.append((value(plain, mask) - sg * value(reverse, mask)) % self.p)
            f, d = np.array(f), np.array(d)
            if check:
                assert np.array_equal(f, d)
            return f, d

        out = 0
        for a, b, phase in subsets(cuts):
            f, d = words(left, a)
            g, e = words(right, b)
            ff, dd = f.reshape(len(f), -1), d.reshape(len(d), -1)
            gg, ee = g.reshape(len(g), -1), e.reshape(len(e), -1)
            gt = g.transpose(0, 2, 1).reshape(len(g), -1).T
            et = e.transpose(0, 2, 1).reshape(len(e), -1).T
            pairs = [dd @ gt - ff @ et]
            for name in ("ad", "jo"):
                gram = self.frames[name][2]
                x = self.coordinates(dd, name) @ gram % self.p
                y = self.coordinates(ff, name) @ gram % self.p
                pairs.append(
                    x @ self.coordinates(gg, name).T - y @ self.coordinates(ee, name).T
                )
            out = out + phase * (np.stack(pairs) % self.p)
            value.cache_clear()
        return np.asarray(out % self.p, dtype=np.uint16).ravel()


@cache
def columns(case, width):
    n = FAMILIES[case][1]
    return np.sort(
        np.random.default_rng(916738 + case).choice(3 * n * n, width, replace=False)
    )


def point(args):
    number, p, cases, width = args
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        e = Evaluation(14337 + number, p)
        return np.concatenate([e.block(case)[columns(case, width)] for case in cases])


def minor(vectors, prime):
    # Exact block Gaussian elimination. Each panel is put in reduced form in
    # its pivot columns so that a single matrix product eliminates it below.
    data = vectors.astype(float)
    rows, cols, det = [], [], 1
    for k in range(0, len(data), 32):
        end = min(k + 32, len(data))
        panel = data[k:end].astype(np.int64)
        active, piv = [], []
        for i in range(len(panel)):
            row = panel[i]
            nonzero = np.flatnonzero(row)
            if not len(nonzero):
                continue
            j = int(nonzero[0])
            lead = int(row[j])
            det = (det * lead) % prime
            unit = row * pow(lead, -1, prime) % prime
            panel = (panel - panel[:, j, None] * unit) % prime
            panel[i] = unit
            active.append(i)
            piv.append(j)
        if active:
            assert np.array_equal(panel[np.ix_(active, piv)], np.eye(len(active)))
            lower = data[end:]
            lower[:] = (lower - lower[:, piv] @ panel[active].astype(float)) % prime
        rows.extend(k + i for i in active)
        cols.extend(piv)
    return rows, cols, det


def gather(numbers, p, cases, width, workers):
    jobs = [(int(i), p, tuple(cases), width) for i in numbers]
    if workers == 1:
        return np.array(list(map(point, jobs)))
    colour(p)
    with ProcessPoolExecutor(workers, mp_context=get_context("fork")) as pool:
        return np.array(list(pool.map(point, jobs, chunksize=2)))


def write(path, value):
    tmp = path.with_suffix(".part")
    tmp.write_text(json.dumps(value, separators=(",", ":")) + "\n")
    tmp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=("search", "replay", "controls"))
    ap.add_argument("--cases", type=int, nargs="+", default=[0, 3, 4])
    ap.add_argument("--points", type=int, default=384)
    ap.add_argument("--width", type=int, default=640)
    ap.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--certificate", type=Path)
    args = ap.parse_args()
    assert __debug__ and args.workers >= 1
    p = args.prime
    if args.action == "controls":
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            e = Evaluation(14337, p, brackets=True)
            for c in args.cases:
                assert not e.block(c, check=True).any()
        print("All selected boundaries vanish on the classical split.")
        return
    if args.action == "search":
        data = gather(range(args.points), p, args.cases, args.width, args.workers)
        rows, cols, det = minor(data, p)
        record = dict(
            prime=p,
            points=args.points,
            cases=args.cases,
            width=args.width,
            rows=rows,
            columns=cols,
            determinant=det,
            rank_lower_bound=len(rows),
        )
        write(args.certificate, record)
    else:
        record = json.loads(args.certificate.read_bytes())
        data = gather(
            record["rows"], p, record["cases"], record["width"], args.workers
        )[:, record["columns"]]
        rows, cols, det = minor(data, p)
        if len(cols) != len(data):
            det = 0
        else:
            det = (
                det
                * (-1) ** sum(x > y for i, x in enumerate(cols) for y in cols[i + 1 :])
                % p
            )
        if p == record["prime"]:
            assert det == record["determinant"]
    print(json.dumps(dict(prime=p, rank=len(rows), determinant=det)), flush=True)


if __name__ == "__main__":
    main()
