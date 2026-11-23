"""Evaluate each return independently by acting on four quadratic defining traces."""
from collections import Counter
from functools import cache
from itertools import product
from pathlib import Path
import json
from jet_algebra import E, K, action, parity, transpose
from verify import WEIGHTS

HERE = Path(__file__).resolve().parent


def odd(letter):
    return (1+letter[2].bit_count()) % 2


@cache
def canonical(trace_factors):
    signed = 1
    factors = []
    for x, y in trace_factors:
        if odd(x) and odd(y):
            if x == y:
                return 0, ()
            if x > y:
                signed = -signed
        factors.append((min(x, y), max(x, y)))
    odd_factors = [t for t in factors if odd(t[0]) != odd(t[1])]
    if len(set(odd_factors)) != len(odd_factors):
        return 0, ()
    signed *= (-1)**sum(t>u for n, t in enumerate(odd_factors) for u in odd_factors[n+1:])
    return signed, tuple(sorted(factors))


@cache
def one(word, op):
    row = Counter()
    prefix = 1
    for i, t in enumerate(word):
        for j, l in enumerate(t):
            n, x = action(l, op)
            if n:
                changed = t[:j] + (x,) + t[j+1:]
                s, w = canonical(word[:i] + (changed,) + word[i+1:])
                if s:
                    row[w] += n*s*prefix
            if odd(l) and parity(op):
                prefix = -prefix
    return tuple((w, n) for w, n in row.items() if n)


def act(v, op):
    row = Counter()
    for w, n in v.items():
        for u, k in one(w, op):
            row[u] += n*k
    return {w: n for w, n in row.items() if n}


def apply(word, v):
    for op in reversed(word):
        v = act(v, op)
    return v


def run():
    pair = ((0, 1, 0), (1, 0, 0))
    terminal = {(pair,) * 4: 1}
    reverse = tuple(E(f,d) for d in range(2) for f in range(2,5))
    primary = apply(reverse, terminal)
    assert primary
    expected = WEIGHTS[0]
    for i,j in product(range(5),repeat=2):
        if i < 2 <= j:
            continue
        v = act(primary, E(i,j))
        factor = expected[i] if i == j else 0
        assert v == ({w:factor*n for w,n in primary.items()} if factor else {})
    assert all(not act(primary, K(i)) for i in range(5))
    back = apply(tuple(transpose(o) for o in reversed(reverse)), primary)
    top_return, = back.values()
    assert top_return and back == {next(iter(terminal)):top_return}
    results = []
    reference = json.loads((HERE/'returned.json').read_bytes())
    for level in reference['levels']:
        count = 0
        for block in level['blocks']:
            words = [tuple(map(tuple, w)) for w in block['words']]
            forward = [apply(w, primary) for w in words]
            for w, values in zip(words, block['matrix']):
                lowering = tuple(transpose(o) for o in reversed(w))
                for v, k in zip(forward, values):
                    assert apply(lowering, v) == ({u:k*n for u,n in primary.items()} if k else {})
                    count += 1
            one.cache_clear()
            canonical.cache_clear()
        results.append(dict(level=level['level'],exact_word_pairings=count))
    return dict(primary_trace_terms=len(primary),sixfold_return=top_return,levels=results)


if __name__ == '__main__':
    assert __debug__
    print(json.dumps(run(), indent=2))
