"""Rational Pfaffian products paired with the native orthogonal one-loop map."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial, lcm
from multiprocessing import get_context
from pathlib import Path
import sys

from bmn_native_pairs import PRIMES, decode, terminal_key


def make_source(line, active, seed):
    _, *values = map(int, line.split())
    matrices = [(a,b,c,-a) for a,b,c in (values[i:i+3] for i in (0,3,6))]
    u, v = values[9:12], values[12:15]

    @lru_cache(None)
    def auxiliary(w):
        a, b, c, d = 1, 0, 0, 1
        for k in w:
            e, f, g, h = matrices[int(k)]
            a, b, c, d = a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h
        return a+d

    @lru_cache(None)
    def logarithm(w):
        raw = {}
        for word, sign in ((w, 1), (w[::-1], (-1)**len(w))):
            for i in range(len(w)):
                t = word[i:]+word[:i]
                assert t not in raw or raw[t] == sign
                raw[t] = sign
        return Fraction((-1)**(len(w)+1)*sum(s*auxiliary(t) for t,s in raw.items()), 2*len(w))

    @lru_cache(None)
    def exponential(factors):
        value = Fraction(1)
        for t, n in Counter(factors).items():
            value *= logarithm(t)**n/factorial(n)
        return value

    def coefficient(w):
        traces = w.split('/')
        t, = [r for r in traces if any(l in '567' for l in r)]
        traces.remove(t)
        value = Fraction()
        for b in set(traces):
            if len(b) == 2:
                i, j = map(int, b)
                bilinear = u[i]*v[j] + (u[j]*v[i] if i != j else 0)
                tail = traces.copy()
                tail.remove(b)
                value += bilinear*exponential(tuple(tail))
        return value*active[t][seed]/(1+sum(len(r)%2 for r in traces))
    return coefficient


def residue(v, p):
    return v.numerator*pow(v.denominator, -1, p) % p


def score_chunk(batch):
    from rank_four_lowered_replay import Q1d_multi
    scores = [Fraction(), Fraction()]
    hits = 0
    for w, checks in batch:
        native = 0
        for raw, n in Q1d_multi(decode(w), 13, 'SO_odd').items():
            sign, key = terminal_key(raw)
            if sign and key in DUAL:
                assert n.denominator == 1
                native += int(n)*sign*DUAL[key]
        coeff = Fraction(native, DEN)
        sources = [fn(w) for fn in SOURCES]
        hits += bool(native and any(sources))
        for p, expected in zip(PRIMES, checks):
            triple = [residue(e, p) for e in [coeff, *sources]]
            assert triple == expected, (p, w, triple, expected)
        for i, v in enumerate(sources):
            scores[i] += coeff*v
    return len(batch), scores, hits


def run(a):
    global DUAL, DEN, SOURCES
    sys.path.insert(0, str(a.repo.resolve()/'algorithms/dual_singlets_20261118'))
    from mixed_incoming import so_canon
    values = {}
    for line in a.dual.read_text().splitlines():
        k, word = line.split()
        sign, w = so_canon(decode(word))
        assert sign and w not in values
        values[w] = sign*Fraction(k)
    DEN = lcm(*(v.denominator for v in values.values()))
    DUAL = {w: int(k*DEN) for w, k in values.items()}
    active = {}
    for line in a.active.read_text().splitlines():
        word, denominator, first, second = line.split()
        active[word] = [Fraction(int(k), int(denominator)) for k in (first, second)]
    lines = a.parameters.read_text().splitlines()
    SOURCES = [make_source(lines[i], active, i) for i in (0, 1)]
    inputs = []
    for path in a.parents:
        entries = {}
        for line in path.read_text().splitlines():
            w, *v = line.split()
            assert w not in entries
            entries[w] = list(map(int, v))
        inputs.append(entries)
    assert inputs[0].keys() == inputs[1].keys()
    rows = [(w, (inputs[0][w],inputs[1][w])) for w in sorted(inputs[0])]
    batches = (rows[i:i+64] for i in range(0, len(rows), 64))
    total, count, hits = [Fraction(), Fraction()], 0, 0
    with ProcessPoolExecutor(a.workers, mp_context=get_context('fork')) as pool:
        for n, v, nz in pool.map(score_chunk, batches):
            count += n
            hits += nz
            for i in range(2):
                total[i] += v[i]
            if count % 8192 == 0:
                print('exact parents', count, flush=True)
    surviving = 358475625*total[0]-832*total[1]
    report = dict(parents=count, contributing=hits, dual_denominator=DEN, active_traces=len(active),
                  dual_words=len(DUAL), corrected_second_on_first=list(map(str,total)),
                  corrected_second_on_combination=str(surviving), third_pairing=str(-surviving),
                  third_pairing_residues=[residue(-surviving,p) for p in PRIMES],
                  ordered_parent_keys_sha256=sha256(''.join(w+'\n' for w,_ in rows).encode()).hexdigest())
    print(report, flush=True)
    a.out.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('repo', 'parameters', 'active', 'dual', 'out'):
        p.add_argument('--'+k, type=Path, required=True)
    p.add_argument('--parents', type=Path, nargs=2, required=True)
    p.add_argument('--workers', type=int, default=2)
    assert __debug__
    run(p.parse_args())
