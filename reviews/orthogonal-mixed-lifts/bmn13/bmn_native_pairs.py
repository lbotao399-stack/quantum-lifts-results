"""Replay every parent pairing with the pinned native Q1d_multi at both primes."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from multiprocessing import get_context
from pathlib import Path
import sys

PRIMES = (1000003, 1000033)
LETTERS = tuple((0, 0, (i,)) for i in range(3)) + ((1, 0, ()), (0, 1, ())) + tuple(
    (0, 0, pair) for pair in ((0, 1), (0, 2), (1, 2)))


def decode(word):
    return tuple(tuple(LETTERS[int(i)] for i in trace) for trace in word.split('/'))


def load_duals(paths):
    from mixed_incoming import so_canon
    dual = {}
    for i, path in enumerate(paths):
        for line in path.read_text().splitlines():
            coefficient, word = line.split()
            sign, key = so_canon(decode(word))
            assert sign and not dual.get(key, [0, 0])[i]
            dual.setdefault(key, [0, 0])[i] = sign * int(coefficient)
    return dual


@lru_cache(16384)
def terminal_key(trace):
    from mixed_incoming import so_canon
    return so_canon(trace)


def native_chunk(rows):
    from rank_four_lowered_replay import Q1d_multi
    scores = [[0, 0], [0, 0]]
    nonzero = [0, 0]
    for word, entries in rows:
        coefficients = [0, 0]
        for trace, value in Q1d_multi(decode(word), 13, 'SO_odd').items():
            sign, key = terminal_key(trace)
            if sign and key in DUAL:
                assert Fraction(value).denominator == 1
                for i, v in enumerate(DUAL[key]):
                    coefficients[i] += int(value) * sign * v
        for i, (p, (expected, first, second)) in enumerate(zip(PRIMES, entries)):
            actual = coefficients[i] % p
            assert actual == expected, (p, word, actual, expected)
            scores[i][0] += first * actual
            scores[i][1] += second * actual
            nonzero[i] += bool(actual and (first or second))
    return len(rows), scores, nonzero


def run(args):
    global DUAL
    sys.path.insert(0, str(args.repo.resolve() / 'algorithms/dual_singlets_20261118'))
    DUAL = load_duals(args.duals)
    inputs = []
    for path in args.parents:
        entries = {}
        for line in path.read_text().splitlines():
            w, *v = line.split()
            assert w not in entries
            entries[w] = tuple(map(int, v))
        inputs.append(entries)
    assert inputs[0].keys() == inputs[1].keys()
    rows = [(w, (inputs[0][w], inputs[1][w])) for w in sorted(inputs[0])]
    if args.limit:
        rows = rows[::max(1, len(rows)//args.limit)][:args.limit]
    scores, nonzero, count = [[0, 0], [0, 0]], [0, 0], 0
    batches = (rows[i:i+64] for i in range(0, len(rows), 64))
    with ProcessPoolExecutor(args.workers, mp_context=get_context('fork')) as pool:
        for n, values, hits in pool.map(native_chunk, batches):
            count += n
            for i, p in enumerate(PRIMES):
                nonzero[i] += hits[i]
                for j in range(2):
                    scores[i][j] = (scores[i][j]+values[i][j]) % p
            if count % 4096 == 0:
                print('native parents', count, flush=True)
    result = dict(parents=count, primes=PRIMES, dual_words=len(DUAL), nonzero=nonzero,
                  source_pairings=scores, parent_keys_sha256=sha256(
                      ''.join(w+'\n' for w, _ in rows).encode()).hexdigest())
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(result, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--duals', type=Path, nargs=2, required=True)
    parser.add_argument('--parents', type=Path, nargs=2, required=True)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert __debug__ and args.workers >= 1
    run(args)
