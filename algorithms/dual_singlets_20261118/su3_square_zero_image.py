"""Restrict the archived first supercharge images to triples with a common image line."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

ALGORITHMS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ALGORITHMS / "su3_j34_20261118"))
from replay_formal_blocks import ALPH, INPUT, LETTER, stream_rows

SCALARS = {LETTER[0, 0, (i,)] for i in range(3)}
ODD = tuple(len(theta) % 2 == 0 for _, _, theta in ALPH)
UNIT = {(b"", 0): 1}


def multiply(left, right):
    """A monomial has commuting coordinates in sorted order and an exterior bit set."""
    out = Counter()
    for (e, mask), a in left.items():
        for (f, other), b in right.items():
            if mask & other:
                continue
            # Each odd variable on the right crosses the larger variables on the left.
            parity = 0
            bits = other
            while bits:
                low = bits & -bits
                parity ^= (mask >> low.bit_length()).bit_count() & 1
                bits -= low
            out[bytes(sorted(e + f)), mask | other] += (-1 if parity else 1) * a * b
    return {m: v for m, v in out.items() if v}


@lru_cache(None)
def entry(k, i, j):
    if k in SCALARS and (i != 0 or j == 0):
        return {}
    position = 3 * i + j
    terms = ((position, 1),) if position < 8 else ((0, -1), (4, -1))
    out = {}
    for pos, sign in terms:
        var = 8*k + pos
        monomial = (b"", 1 << var) if ODD[k] else (bytes([var]), 0)
        out[monomial] = sign
    return out


def adjacent_scalars(trace):
    return any(k in SCALARS and trace[i-1] in SCALARS for i, k in enumerate(trace))


@lru_cache(maxsize=18000)
def trace_polynomial(trace):
    if adjacent_scalars(trace):
        return {}
    first_scalar = next((i for i, k in enumerate(trace) if k in SCALARS), None)
    sign = 1
    starts = range(3)
    if first_scalar is not None:
        left, right = trace[:first_scalar], trace[first_scalar:]
        sign = (-1)**(sum(ODD[k] for k in left) * sum(ODD[k] for k in right))
        trace = right + left
        starts = (0,)
    result = Counter()
    for start in starts:
        paths = {start: UNIT}
        for place, k in enumerate(trace):
            next_paths = defaultdict(Counter)
            for i, previous in paths.items():
                for j in ((start,) if place == len(trace)-1 else range(3)):
                    next_paths[j].update(multiply(previous, entry(k, i, j)))
            paths = next_paths
        result.update(paths.get(start, {}))
    return {m: sign*v for m, v in result.items() if v}


def coordinate(var):
    return [var // 8, var % 8]


def describe(key, value):
    even, mask = key
    odd = [i for i in range(mask.bit_length()) if (mask >> i) & 1]
    return dict(even=list(map(coordinate, even)), odd=list(map(coordinate, odd)), coefficient=str(value))


def expand_content(rows):
    coefficients = Counter()
    for word, scale in rows:
        polynomial = UNIT
        for trace in word:
            polynomial = multiply(polynomial, trace_polynomial(trace))
            if not polynomial:
                break
        for monomial, multiplicity in polynomial.items():
            coefficients[monomial] += scale * multiplicity
    return {m: v for m, v in coefficients.items() if v}


def centralizers(terms, denominator):
    gauginos = {LETTER[1, 0, ()], LETTER[0, 1, ()]}
    rank_two = Counter()
    rank_one = Counter()
    for (even, mask), n in terms:
        indices = [i for i in range(mask.bit_length()) if mask & (1 << i)]
        slots = [i % 8 for i in indices if i // 8 in gauginos]
        if all(s in (1, 2) for s in slots):
            rank_two[even, mask] += n
        if any(i // 8 in SCALARS and i % 8 != 1 for i in even) or any(s in (3, 5, 6) for s in slots):
            continue
        reduced = [i-4 if i // 8 in gauginos and i % 8 == 4 else i for i in indices]
        if len(set(reduced)) == len(reduced):
            sign = (-1)**sum(i > j for p,i in enumerate(reduced) for j in reduced[p+1:])
            rank_one[even, sum(1 << i for i in reduced)] += sign*n
    result = {}
    for name, sparse in (("two_independent_scalars", rank_two), ("one_independent_scalar", rank_one)):
        nonzero = sorted((m, n) for m,n in sparse.items() if n)
        result[name] = dict(nonzero_monomials=len(nonzero),
                            first_coefficient=describe(nonzero[0][0], Fraction(nonzero[0][1], denominator)) if nonzero else None)
        if name == "one_independent_scalar":
            result[name]["polynomial"] = [describe(m, Fraction(n, denominator)) for m,n in nonzero]
    return result


def calculate(level, input_path=None, workers=1, reduce_gauginos=False):
    if level == 32:
        folder = Path(__file__).resolve().parent / "taylor_hooks/su3_level32"
        path = folder / "Q1.jsonl.gz"
        record = json.loads((folder / "order1.json").read_bytes())
        expected, total, denominator, orientation = (record["output_sha256"], record["direct_terms"],
                                                     record["common_denominator"], 1)
    else:
        path = ALGORITHMS / "su3_j34_20261118/negative_Q1.jsonl.gz"
        expected, total, denominator, orientation = INPUT, 517468, 120, -1
    surviving_words = Counter()
    contents = defaultdict(list)
    digest = hashlib.sha256()
    count = 0
    for word, value in stream_rows(input_path or path, digest):
        count += 1
        if any(adjacent_scalars(trace) for trace in word):
            continue
        nscalar = sum(k in SCALARS for trace in word for k in trace)
        surviving_words[nscalar] += 1
        scale = Fraction(value) * denominator * orientation
        if scale.denominator != 1:
            raise ValueError("first image exceeds its archived common denominator")
        contents[bytes(sorted(b"".join(word)))].append((word, scale.numerator))
    if count != total or digest.hexdigest() != expected:
        raise ValueError("first image differs from its archive")
    groups = sorted(contents.values(), key=len, reverse=True)
    if workers == 1:
        terms = sorted(item for group in groups for item in expand_content(group).items())
    else:
        with ProcessPoolExecutor(workers) as pool:
            terms = sorted(item for result in pool.map(expand_content, groups) for item in result.items())
    polynomial_hash = hashlib.sha256()
    by_degree = Counter()
    for key, n in terms:
        record = describe(key, Fraction(n, denominator))
        by_degree[sum(v[0] in SCALARS for v in record["even"])] += 1
        polynomial_hash.update((json.dumps(record, separators=(",", ":"))+"\n").encode())
    report = dict(level=level, direct_image="Q1", input_sha256=expected, input_terms=count,
                surviving_cyclic_words_by_scalar_count=dict(sorted(surviving_words.items())),
                terms_by_scalar_count=dict(sorted(by_degree.items())), nonzero_monomials=len(terms),
                expanded_polynomial_sha256=polynomial_hash.hexdigest(),
                first_coefficient=describe(terms[0][0], Fraction(terms[0][1], denominator)) if terms else None)
    if reduce_gauginos:
        report["gauginos_in_scalar_centralizer"] = centralizers(terms, denominator)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, choices=(32, 34), required=True)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--centralizers", action="store_true")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    data = calculate(args.level, args.input, args.workers, args.centralizers)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print(args.level, data["surviving_cyclic_words_by_scalar_count"], data["nonzero_monomials"], flush=True)
