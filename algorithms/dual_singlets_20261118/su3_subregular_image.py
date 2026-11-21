"""Expand the signed SU(3) first image on the subregular scalar chart over Q."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path


ALGORITHMS = Path(__file__).resolve().parents[1]
ARCHIVES = {
    32: ("dual_singlets_20261118/taylor_hooks/su3_level32/Q1.jsonl.gz", 331201, 3,
         "f7810267d100f7cc9c9b4c93e045540ab2291df8b8bb35eeb4ad77feddd47a8a"),
    34: ("su3_j34_20261118/negative_Q1.jsonl.gz", 517468, 120,
         "9435489d5e21f083ac22cb9bfd92b0dd832192252ee8b0a3ac04f656f5c85f14"),
}
FIELDS = "pquradfg"
THETA_DEGREE = (0, 1, 1, 2)
ODD = tuple(bool((k // 16 >= 4) ^ (THETA_DEGREE[k % 4] % 2)) for k in range(128))
UNIT = {b"": 1}
INFINITY = 100


def multiply(left, right, limit):
    result = Counter()
    for x, c in left.items():
        for y, d in right.items():
            if len(x) + len(y) > limit:
                continue
            phase, monomial = ordered(x, y)
            if phase:
                result[monomial] += phase*c*d
    return {m: v for m, v in result.items() if v}


@lru_cache(60000)
def ordered(x, y):
    sign = 1
    for i in x:
        if ODD[i]:
            for j in y:
                if i == j:
                    return 0, b""
                if i > j and ODD[j]:
                    sign = -sign
    return sign, bytes(sorted(x+y))


def generator(field, index):
    if field in "adfg" and index == 0:
        return {}
    return {bytes([16*FIELDS.index(field) + index]): 1}


@lru_cache(None)
def coefficient(factors, index):
    if len(factors) == 1:
        return generator(factors, index)
    if not factors:
        return UNIT if index == 0 else {}
    left, right = factors[0], factors[1:]
    result = Counter()
    # Each Taylor index packs the two spacetime powers and the two remaining odd slots.
    for i in range(16):
        if i & ~index:
            continue
        j = index ^ i
        theta_i, theta_j = i & 3, j & 3
        parity = (left in "adfg") ^ (THETA_DEGREE[theta_i] % 2)
        crossings = bool(theta_i & 2 and theta_j & 1)
        sign = (-1)**(crossings + parity*THETA_DEGREE[theta_j])
        for monomial, n in multiply(generator(left, i), coefficient(right, j), len(factors)).items():
            result[monomial] += sign*n
    return {m: v for m, v in result.items() if v}


B_ENTRIES = ((), ((1, "p"),), ((1, "u"),),
             ((1, ""),), ((1, "q"),), (),
             (), ((1, "r"),), ((-1, "q"),))
A_ENTRIES = (((1, "a"),), ((1, "dp"), (1, "gu")), ((1, "du"), (-2, "fq")),
             ((1, "d"),), ((1, "a"), (1, "dq")), ((1, "f"),),
             ((1, "g"),), ((1, "dr"), (-1, "gq")), ((-2, "a"), (-1, "dq")))


@lru_cache(None)
def entries(letter, vanish=""):
    index = letter >> 1
    result = []
    for terms in B_ENTRIES if letter & 1 else A_ENTRIES:
        polynomial = Counter()
        for scale, factors in terms:
            if any(field in vanish for field in factors):
                continue
            for key, n in coefficient(factors, index).items():
                polynomial[key] += scale*n
        result.append({m: n for m, n in polynomial.items() if n})
    return result


@lru_cache(None)
def costs(letter, vanish=""):
    return [min(map(len, poly), default=INFINITY) for poly in entries(letter, vanish)]


@lru_cache(100000)
def trace_bound(trace, vanish=""):
    # A closed colour path is a summand of the trace. Cancellations can only raise its degree.
    distances = list(costs(trace[0], vanish))
    for letter in trace[1:]:
        next_costs = costs(letter, vanish)
        distances = [min(distances[3*i+k] + next_costs[3*k+j] for k in range(3))
                     for i in range(3) for j in range(3)]
    return min(distances[0], distances[4], distances[8])


@lru_cache(12000)
def trace_polynomial(trace, limit, vanish=""):
    result = Counter()
    for first in range(3):
        suffix = [[0 if k == first else INFINITY for k in range(3)]]
        for letter in reversed(trace):
            previous = suffix[-1]
            cost = costs(letter, vanish)
            suffix.append([min(cost[3*k+j]+previous[j] for j in range(3)) for k in range(3)])
        if suffix[-1][first] > limit:
            continue
        paths = {first: UNIT}
        for i, letter in enumerate(trace):
            following = defaultdict(Counter)
            for k, poly in paths.items():
                for j in (range(3) if i+1 < len(trace) else (first,)):
                    entry = entries(letter, vanish)[3*k+j]
                    cap = limit - suffix[len(trace)-i-1][j]
                    if entry and cap >= costs(letter, vanish)[3*k+j]:
                        following[j].update(multiply(poly, entry, cap))
            paths = {j: poly for j, poly in following.items() if poly}
        result.update(paths.get(first, {}))
    return {m: n for m, n in result.items() if n}


def source(level, archive, limit, vanish=""):
    relative, size, denominator, sha = ARCHIVES[level]
    digest = hashlib.sha256()
    counts = Counter()
    selected = []
    with gzip.open(archive or ALGORITHMS/relative, "rb") as stream:
        for line in stream:
            digest.update(line)
            value, raw = json.loads(line)
            word = tuple(bytes((2*z+w)*8 + sum(1 << i for i in theta) for z, w, theta in tr) for tr in raw)
            if sum(map(len, word)) != 11 or sum(k & 1 for trace in word for k in trace) != 7:
                raise ValueError("expected seven B coefficients and four A coefficients")
            # Only the undifferentiated first scalar has a constant component.
            if 11 - sum(t.count(1) for t in word) > limit:
                counts["scalar_exclusion"] += 1
                continue
            lower = sum(trace_bound(t, vanish) for t in word)
            if lower > limit:
                counts["colour_exclusion"] += 1
                continue
            numerator = Fraction(value) * denominator
            assert numerator.denominator == 1
            counts[str(lower)] += 1
            selected.append([int(numerator), [t.hex() for t in word]])
    if digest.hexdigest() != sha or sum(counts.values()) != size:
        raise ValueError("first image does not match the archive")
    return dict(level=level, sign=1 if level==32 else -1, denominator=denominator, degree=limit,
                full_image_sha256=sha, full_image_words=size, bounds=dict(sorted(counts.items())),
                selected=selected, **({"vanishing_fields":vanish} if vanish else {}))


def expand_words(selected, limit, vanish=""):
    total = Counter()
    for numerator, encoded in selected:
        word = [bytes.fromhex(t) for t in encoded]
        bounds = [trace_bound(t, vanish) for t in word]
        lower = remaining = sum(bounds)
        polynomial = UNIT
        for t, bound in zip(word, bounds):
            remaining -= bound
            tr = trace_polynomial(t, limit-lower+bound, vanish)
            polynomial = multiply(polynomial, tr, limit-remaining)
            if not polynomial:
                break
        for m, n in polynomial.items():
            total[m] += numerator*n
    return {m: n for m, n in total.items() if n}


def expand(packet, workers=1):
    selected, limit = packet["selected"], packet["degree"]
    vanish = packet.get("vanishing_fields","")
    if workers == 1:
        total = expand_words(selected, limit, vanish)
    else:
        total = Counter()
        size = 4000
        with ProcessPoolExecutor(workers) as pool:
            pieces = [pool.submit(expand_words, selected[i:i+size], limit, vanish) for i in range(0, len(selected), size)]
            for piece in pieces:
                total.update(piece.result())
    rows = [[m.hex(), n] for m, n in sorted(total.items()) if n]
    degrees = Counter(len(m)//2 for m, _ in rows)
    return {**{k: v for k, v in packet.items() if k != "selected"}, "selected_words":len(packet["selected"]),
            "nonzero_by_degree":{str(k):degrees[k] for k in range(limit+1)}, "coefficients":rows}


def write(path, data):
    raw = (json.dumps(data, separators=(",", ":"))+"\n").encode()
    path.write_bytes(gzip.compress(raw, mtime=0) if path.suffix == ".gz" else raw)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("prepare", "expand"))
    parser.add_argument("--level", type=int, choices=tuple(ARCHIVES))
    parser.add_argument("--degree", type=int)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires assertions")
    if args.task == "prepare":
        if args.level is None or args.degree is None or args.degree < 0:
            parser.error("prepare requires --level and a nonnegative --degree")
        data = source(args.level, args.input, args.degree)
        print(data["bounds"], flush=True)
    else:
        if args.input is None:
            parser.error("expand requires a compressed packet")
        if args.workers < 1:
            parser.error("the worker count must be positive")
        data = expand(json.loads(gzip.decompress(args.input.read_bytes())), args.workers)
        print(data["nonzero_by_degree"], flush=True)
    write(args.out, data)
