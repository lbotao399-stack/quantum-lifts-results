"""Polarize products of two finite F4 traces in the (4,4,3,3) classical boundary test."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from functools import cache
import json
from pathlib import Path

import numpy as np

from flag_boundaries import HERE, adjoint, flags, polynomial
from single_trace_flag import catalogue, extend, native26, seed_matrix, short_values, trace_values


def necklace(words, length):
    mask = (1 << (2*length))-1
    rev = np.zeros_like(words)
    for j in range(length):
        rev = (rev << 2) + ((words >> (2*j)) & 3)

    def least_rotation(raw):
        current = raw
        best = raw.copy()
        for _ in range(1, length):
            current = ((current << 2) & mask) + (current >> (2*length-2))
            best = np.minimum(best, current)
        return best

    forward, backward = least_rotation(words), least_rotation(rev)
    sign = np.sign(backward-forward) if length % 2 else 1
    return np.minimum(forward, backward), sign


def double_terms(length, seed):
    slots, signs = flags()
    perm = np.random.Generator(np.random.PCG64(143429+1031*seed+length)).permutation(14)

    def word(labels):
        value = slots[labels[0]]
        for k in labels[1:]:
            value = 4*value + slots[k]
        return necklace(value, len(labels))

    long, left_sign = word(perm[:length])
    short, right_sign = word(perm[length:])
    packed = (long.astype(np.int64) << 14) + short
    keys, inverse = np.unique(packed, return_inverse=True)
    coeff = np.bincount(inverse, weights=signs*left_sign*right_sign).astype(np.int64)
    return keys[coeff != 0], coeff[coeff != 0]


class DoubleTraces:
    def __init__(self, length, seeds):
        self.length = length
        terms = [double_terms(length, seed) for seed in seeds]
        keys = np.unique(np.concatenate([k for k, _ in terms]))
        self.indices = np.concatenate([np.searchsorted(keys, k) for k, _ in terms])
        ends = np.cumsum([0] + [len(k) for k, _ in terms])
        self.active = np.flatnonzero(np.diff(ends))
        self.starts = ends[self.active]
        self.coefficients = np.concatenate([v for _,v in terms])
        self.size = len(seeds)
        self.long, self.left = np.unique(keys >> 14, return_inverse=True)
        self.short, self.right = np.unique(keys & ((1<<14)-1), return_inverse=True)

    def entries(self, long, short, prime):
        values = long[self.left].astype(np.int64) * short[self.right]
        out = np.zeros(self.size, dtype=np.int64)
        if len(self.active):
            terms = self.coefficients*(values[self.indices] % prime)
            out[self.active] = np.add.reduceat(terms,self.starts) % prime
        return out


def powers(x, prime):
    dim = x.shape[1]
    table = {1: x}
    for n in range(2, 8):
        left = table[n-1].reshape(-1, dim)
        value = np.empty((4**n, dim, dim), dtype=np.float64)
        for j in range(4):
            # Accelerate can leave stale floating-point flags on exact GEMMs.
            # Each dot product sums dim nonnegative integer terms.
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                part = left @ x[j]
            assert 0 <= part.min() and part.max() <= dim*(prime-1)**2
            part -= np.floor(part/prime)*prime
            value[j::4] = part.reshape(-1, dim, dim)
        table[n] = value
    return table


@cache
def family(name):
    if name == "26":
        return native26()
    if name == "52":
        return adjoint()
    raise ValueError(name)


def default_rows(count):
    previous = json.loads((HERE/"flag_433.json").read_bytes())["rows"]
    extra = []
    for length in range(7, 13):
        for a in ("26", "52"):
            for b in (("26",) if length >= 11 else ("26", "52")):
                if length == 7 and a == "52" and b == "26":
                    continue
                extra.extend([f"{length}.{a}.{b}", seed] for seed in range(count))
    return previous + extra


class Evaluation:
    def __init__(self, rows):
        self.rows = rows
        singles = sorted({s for name, s in rows if name in ("26", "52")})
        shorts = [s for name, s in rows if name == "short"]
        self.single = seed_matrix(singles) if singles else None
        self.short_terms = [polynomial(s) for s in shorts]
        index = {key: dict(zip(seq, range(len(seq)))) for key, seq in
                 (("26", singles), ("52", singles), ("short", shorts))}
        self.doubles = {}
        self.pairs = {}
        for name, _ in rows:
            if "." not in name or name in self.pairs:
                continue
            n, a, b = name.split(".")
            self.pairs[name] = (int(n), a, b)
        for n in sorted({p[0] for p in self.pairs.values()}):
            seeds = sorted({s for name, s in rows if name in self.pairs and self.pairs[name][0] == n})
            self.doubles[n] = DoubleTraces(n, seeds)
            lookup = dict(zip(seeds, range(len(seeds))))
            for name, (length, _, _) in self.pairs.items():
                if n == length:
                    index[name] = lookup
        self.positions = [index[name][s] for name, s in rows]
        self.names = sorted({k for k, _ in rows})
        for name in ("26", "52"):
            family(name)

    def point(self, number, prime):
        coefficients = np.random.Generator(np.random.PCG64(144336+number)).integers(-3, 4, size=(4, 52))
        data, left, right = {}, {}, {}
        for name in ("26", "52"):
            rep = family(name)
            dim = rep.shape[1]
            assert dim**2*(prime-1)**2 < 2**53
            x = ((coefficients @ rep.reshape(52,-1)).reshape(4,dim,dim) % prime).astype(np.float64)
            table = powers(x, prime)
            if name in self.names:
                entries, columns, starts = self.single
                values = trace_values(table, catalogue(), 14, prime)
                data[name] = np.add.reduceat(entries*values[columns], starts) % prime
            if name == "52" and self.short_terms:
                data["short"] = short_values(table, self.short_terms, prime)
            for n, double in self.doubles.items():
                left[n,name] = trace_values(table, double.long, n, prime)
                right[n,name] = trace_values(table, double.short, 14-n, prime)
            del table
        for name, (n, a, b) in self.pairs.items():
            data[name] = self.doubles[n].entries(left[n,a], right[n,b], prime)
        return np.array([data[name][j] for (name, _), j in zip(self.rows, self.positions)], dtype=np.int64)


def calculate(args):
    ref = json.loads(args.replay.read_bytes()) if args.replay else None
    rows = ref["rows"] if ref else default_rows(args.seeds)
    points = ref["points"] if ref else range(args.points)
    evaluator = Evaluation(rows)
    basis, kept, recipes = {}, [], []
    det = 1

    def result():
        return dict(prime=args.prime, rank=len(kept), determinant=det, points=kept, rows=recipes)

    with ThreadPoolExecutor(args.workers) as pool:
        for lo in range(0, len(points), args.workers):
            batch = points[lo:lo+args.workers]
            for j, values in zip(batch, pool.map(lambda i: evaluator.point(i, args.prime), batch)):
                pivot = extend(basis, values, args.prime)
                if pivot:
                    k, value = pivot
                    kept.append(j)
                    recipes.append(rows[k])
                    det = det*value % args.prime
            if args.checkpoint:
                tmp = args.checkpoint.with_suffix(".part")
                tmp.write_text(json.dumps(result())+"\n")
                tmp.replace(args.checkpoint)
            if lo//args.workers % 4 == 0:
                print("points", lo+len(batch), "rank", len(kept), flush=True)
            assert len(kept) <= 703
            if len(kept) == 703:
                break
    out = result()
    if ref:
        assert kept == ref["points"] and recipes == ref["rows"]
        if ref["prime"] == args.prime:
            assert ref == out
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=120)
    parser.add_argument("--points", type=int, default=704)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--prime", type=int, choices=(1009,1013), default=1009)
    parser.add_argument("--replay", type=Path)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or min(args.points,args.seeds,args.workers) < 1:
        parser.error("positive sizes and enabled assertions are required")
    record = calculate(args)
    args.out.write_text(json.dumps(record, indent=2)+"\n")
    print("rank", record["rank"], "determinant", record["determinant"], flush=True)
