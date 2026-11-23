"""Exact native Pfaffian expansion by graded cycle decomposition.

This is a separate experimental engine.  Coefficients are compared with
path_relation_completion.expand_paths; no physical normalization is changed.
"""
from pathlib import Path
import os
import sys
from collections import Counter
from functools import lru_cache
from itertools import product

CODE_DIR = Path(__file__).resolve().parent
HERE = Path(os.environ['FORTUITY_RUN_DIR']).resolve() if os.environ.get('FORTUITY_RUN_DIR') else CODE_DIR
BASE_PACKAGE=Path(os.environ['FORTUITY_BASE_PACKAGE']).resolve() if os.environ.get('FORTUITY_BASE_PACKAGE') else None
sys.path.insert(0, str(BASE_PACKAGE/'symplectic/sp5' if BASE_PACKAGE else HERE.parent/'fortuity_discovery_20261123/jets'))
import path_relation_completion as pc
import pfaffian_relations as pf


def clean(row):
    return {w: c for w, c in row.items() if c}


def cache_native_canonicalization():
    """Memoize unchanged upstream functions; no convention is altered."""
    if not hasattr(pf.m.fw.canon_word,'cache_info'):
        pf.m.fw.canon_word=lru_cache(maxsize=262144)(pf.m.fw.canon_word)
    if not hasattr(pf.m.so_canon,'cache_info'):
        pf.m.so_canon=lru_cache(maxsize=65536)(pf.m.so_canon)


@lru_cache(maxsize=32768)
def cycle(paths):
    """One anchored, directed cycle, including all path orientations."""
    orientations = []
    for path in paths:
        odd = sum(map(pc.is_ferm, path))
        sign = (-1) ** (len(path) + odd * (odd - 1) // 2)
        row = Counter({path: 1})
        row[path[::-1]] += sign
        orientations.append(tuple(clean(row).items()))
    result = Counter()
    for choice in product(*orientations):
        word = tuple(l for path, _ in choice for l in path)
        coefficient = 1
        for _, c in choice:
            coefficient *= c
        sign, canonical = pf.m.so_canon((word,))
        if sign:
            result[canonical] -= coefficient * sign
    result = clean(result)
    # A cycle and its reversed order share the factor 1/2.  For length >= 3
    # the division is performed only after all anchored orders are summed.
    return result


def multiply(left, right):
    result = Counter()
    for u, c in left.items():
        for v, d in right.items():
            sign, word = pf.m.so_canon(u + v)
            if sign:
                result[word] += c * d * sign
    return clean(result)


def expand_paths(paths, rank=None):
    """Return exactly the original endpoint-antisymmetrizer polynomial."""
    paths = tuple(paths)
    if rank is not None:
        assert len(paths) == rank + 1
    input_sign, paths = pc.canonical_paths(paths)
    if not input_sign:
        return {}

    @lru_cache(None)
    def covers(remaining):
        if not remaining:
            return {(): 1}
        answer = Counter()

        def grow(sequence, rest, coefficient):
            row = multiply(cycle(sequence), covers(rest))
            for w, c in row.items():
                answer[w] += coefficient * c
            counts = Counter(rest)
            odd_before = 0
            for path, count in sorted(counts.items()):
                odd = pc.path_odd(path)
                if odd:
                    assert count == 1
                index = rest.index(path)
                left = rest[:index] + rest[index + 1:]
                sign = (-1) ** (odd * odd_before)
                grow(sequence + (path,), left, coefficient * count * sign)
                odd_before += odd * count

        grow((remaining[0],), remaining[1:], 1)
        answer = clean(answer)
        assert all(c % 2 == 0 for c in answer.values())
        return {w: c // 2 for w, c in answer.items()}

    # Endpoint vertices precede internal vertices in the native graph recipe.
    # Each nontrivial internal chain has exactly one descending fixed edge.
    overall = input_sign * (-1) ** sum(len(path) > 2 for path in paths)
    return {w: overall * c for w, c in covers(paths).items()}


@lru_cache(maxsize=4096)
def graph_components(length, rank, remaining):
    marked = 2 * rank + 2
    adjacent = {a: b for a, b in remaining} | {b: a for a, b in remaining}
    seen = set()
    closing = []
    for start in range(marked):
        if start in seen:
            continue
        seen.add(start)
        half = 2 * start + 1
        while True:
            half = adjacent[half]
            vertex = half // 2
            seen.add(vertex)
            if vertex < marked:
                closing.append(tuple(sorted((2 * start, 2 * vertex))))
                break
            half ^= 1
    closing.sort()
    crossings = sum(a < c < b < d for i, (a,b) in enumerate(closing)
                    for c,d in closing[i+1:])
    graph_sign, traces, order = pf.graph_traces(length, tuple(closing) + remaining)
    paths = tuple(t for t in traces if min(t) < marked)
    cycles = tuple(t for t in traces if min(t) >= marked)
    assert traces == paths + cycles and len(paths) == rank + 1
    canonical_sign = (-1) ** (rank + 1 + sum(len(p) > 2 for p in paths))
    return paths, cycles, order, (-1) ** crossings * graph_sign * canonical_sign


def identity(letters, rank, remaining):
    paths, cycles, order, sign = graph_components(len(letters), rank, remaining)
    fermions = [i for i in order if pc.is_ferm(letters[i])]
    sign *= (-1) ** sum(a > b for j,a in enumerate(fermions) for b in fermions[j+1:])
    row = expand_paths(tuple(tuple(letters[i] for i in p) for p in paths), rank)
    closed = tuple(tuple(letters[i] for i in c) for c in cycles)
    result = Counter()
    for word,c in row.items():
        phase, word = pf.m.so_canon(word + closed)
        if phase:
            result[word] += sign * phase * c
    return clean(result)


if __name__ == '__main__':
    import argparse
    import json
    import pickle
    import random
    import time
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples', type=int, default=20)
    args = parser.parse_args()
    rng = random.Random(73193)
    letters = ((0,0,(0,)), (0,0,(1,)), (1,0,(0,)),
               (1,0,()), (0,1,()), (2,0,()))
    results = []
    for rank in range(4):
        for trial in range(args.samples):
            paths = tuple(tuple(rng.choice(letters) for _ in range(rng.choice((2,2,3,4))))
                          for _ in range(rank + 1))
            t = time.monotonic()
            old = pc.expand_paths(paths, rank)
            t1 = time.monotonic()
            new = expand_paths(paths, rank)
            t2 = time.monotonic()
            if old != new:
                bad = next(w for w in old.keys() | new.keys() if old.get(w,0) != new.get(w,0))
                print(json.dumps(dict(rank=rank,trial=trial,paths=paths,
                                      word=bad,old=old.get(bad,0),new=new.get(bad,0))),flush=True)
                raise AssertionError('Exact coefficient mismatch')
            results.append(dict(rank=rank,trial=trial,terms=len(new),
                                endpoint_seconds=t1-t,cycle_seconds=t2-t1))
        print('verified',rank,args.samples,flush=True)
    (HERE / 'cycle_engine_validation.json').write_text(json.dumps(dict(
        exact_native_coefficients=True,seed=73193,results=results),indent=2)+'\n')
