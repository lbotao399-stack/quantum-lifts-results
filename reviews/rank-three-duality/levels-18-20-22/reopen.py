"""Rebuild each selected Cartan minor over the integers, then take two field determinants."""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
from math import prod
import argparse
import gzip
import hashlib
import json

from cartan_sectors import products, read_census
from graviton import cartan
from verify import WEIGHTS, charges, highest, shift

HERE = Path(__file__).resolve().parent


def odd(letter):
    return (1 + letter[2].bit_count()) % 2


def exterior(items, parities):
    seen = set()
    sign = 1
    for letter, parity in zip(items, parities):
        if parity:
            if letter in seen:
                return 0, ()
            sign *= (-1)**sum(x > letter for x in seen)
            seen.add(letter)
    return sign, tuple(sorted(items))


@lru_cache(maxsize=None)
def one_axis(g):
    """Keep original trace slots labelled until all outer derivatives have acted."""
    letters = tuple((0, 0, 1 << i) for i, n in enumerate(g[2:5]) for _ in range(n))
    letters += ((1, 0, 0),)*g[0] + ((0, 1, 0),)*g[1]
    terms = {letters: 2}
    axes = [j+2 for j in range(3) if g[7] & (1 << j)] + [0]*g[5] + [1]*g[6]
    for axis in axes:
        out = Counter()
        for word, coefficient in terms.items():
            preceding = 0
            for j, letter in enumerate(word):
                p, q, f = letter
                if axis < 2:
                    factor = letter[axis] + 1
                    new = (p + int(axis == 0), q + int(axis == 1), f)
                else:
                    mask = 1 << (axis-2)
                    factor = 0 if mask & f else (-1)**(preceding + (f >> (axis-1)).bit_count())
                    new = (p, q, mask | f)
                if factor:
                    out[word[:j] + (new,) + word[j+1:]] += coefficient * factor
                preceding += odd(letter)
        terms = {word: n for word, n in out.items() if n}
    out = Counter()
    for word, n in terms.items():
        s, key = exterior(word, map(odd, word))
        if s:
            out[key] += s*n
    return tuple((key, v) for key, v in out.items() if v)


@lru_cache(maxsize=4800)
def on_one_colour(factors):
    result = Counter()
    for choices in product(*(one_axis(g) for g in factors)):
        letters = tuple(l for word, _ in choices for l in word)
        s, w = exterior(letters, map(odd, letters))
        if s:
            result[w] += s * prod(v for _, v in choices)
    return tuple((word, n) for word, n in result.items() if n)


@lru_cache(maxsize=None)
def partitions(count):
    assignments = []
    for colours in product(range(3), repeat=count):
        used = tuple(dict.fromkeys(colours))
        if used != tuple(range(len(used))):
            continue
        blocks = tuple(tuple(i for i, c in enumerate(colours) if c == h) for h in used)
        assignments.append((colours, blocks, (1, 3, 6, 6)[len(used)]))
    return assignments


def expand(factors):
    parities = [(g[0]+g[1]+g[7].bit_count()) % 2 for g in factors]
    result = Counter()
    for colours, blocks, multiplicity in partitions(len(factors)):
        inversions = sum(colours[i] > colours[j] for i, a in enumerate(parities) if a
                         for j in range(i+1, len(factors)) if parities[j])
        phase = multiplicity * (-1)**inversions
        groups = [tuple(factors[i] for i in block) for block in blocks]
        gp = [sum(parities[i] for i in block) % 2 for block in blocks]
        for choices in product(*(on_one_colour(g) for g in groups)):
            slots = tuple(w for w, _ in choices) + ((),)*(3-len(choices))
            s, key = exterior(slots, gp + [0]*(3-len(choices)))
            if s:
                result[key] += phase*s*prod(v for _, v in choices)
    return {key: v for key, v in result.items() if v}


def determinant(rows, prime):
    """Fixed column order, signed row swaps and ordinary row elimination."""
    mat = [{j: v % prime for j, v in row.items() if v % prime} for row in rows]
    value = 1
    for j in range(len(mat)):
        i = next((k for k in range(j, len(mat)) if j in mat[k]), None)
        if i is None:
            return 0
        if i != j:
            mat[i], mat[j] = mat[j], mat[i]
            value = -value
        t = mat[j][j]
        value = value * t % prime
        u = pow(t, -1, prime)
        for row in mat[j+1:]:
            if j not in row:
                continue
            multiple = row[j]*u
            for h, x in mat[j].items():
                y = (row.get(h, 0) - multiple*x) % prime
                if y:
                    row[h] = y
                else:
                    row.pop(h, None)
    return value


def matrix(record):
    """The first-seen integer IDs encode coordinates only, never coefficients."""
    column_ids = {}
    selected = set(record['selected_rows'])
    columns = {k: j for j, k in enumerate(record['selected_columns'])}
    rows = []
    weight = tuple(record['weight'])
    for i, factors in enumerate(products(weight, record['length'])):
        original = cartan(factors)
        for key in original:
            column_ids.setdefault(key, len(column_ids))
        if i in selected:
            values = expand(factors)
            assert values == original
            rows.append({columns[column_ids[k]]: v for k, v in values.items() if column_ids[k] in columns})
        if i == max(selected):
            break
    n = len(rows)
    assert n == len(columns) == record['cartan_rank'] == record['limit']
    results = [determinant(rows, p) for p in (1009, 1013)]
    assert results[0] == record['determinant'] and all(results)
    cartan.cache_clear()
    on_one_colour.cache_clear()
    return dict(level=record['level'], length=record['length'], dotted=record['dotted'],
                weight=weight, order=n, determinants=results)


def upper_bounds(group, j):
    raw = gzip.decompress((HERE/f'{group}_full_j{j}.jsonl.gz').read_bytes())
    rows = [json.loads(line) for line in raw.splitlines()]
    ranks = {(tuple(r['sec']), r['Y']): r['r_up'] for r in rows}
    by_sector = Counter()
    for r in rows:
        sec, y = tuple(r['sec']), r['Y']
        predecessor = (sec[0]-1, sec[1]-1, *sec[2:])
        h = r['m'] - ranks[sec, y] - ranks.get((predecessor, y+1), 0)
        assert h >= 0 and r['check']
        if h:
            by_sector[sec, y] = h
    return by_sector


def check_census(records):
    provenance = json.loads((HERE/'provenance.json').read_bytes())
    for path, expected in provenance['blobs'].items():
        name = Path(path).name
        raw = (gzip.decompress((HERE/(name+'.gz')).read_bytes()) if name.endswith('.jsonl')
               else (HERE/name).read_bytes())
        assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest() == expected
    levels = json.loads((HERE/'returned.json').read_bytes())['levels']
    verified = []
    for level in levels:
        j = level['level']
        b, c = (upper_bounds(g, j) for g in ('so7','usp6'))
        delta = b.copy()
        delta.subtract(c)
        differences = {k: v for k, v in delta.items() if v}
        comps = [Counter(), Counter()]
        for block in level['blocks']:
            for word in block['words']:
                w = shift(word)
                for seed, bank in zip(WEIGHTS, comps):
                    bank[charges(tuple(a+b for a,b in zip(seed,w)))] += 1
        returned = Counter()
        for bank in comps:
            returned.update({(q, 8): v for q, v in highest(bank).items()})
        assert differences == returned
        expected = set()
        for y in range(1, j//2+1):
            wb, wc = (read_census(g, j, y) for g in ('so7','usp6'))
            for d in range(j%2, j//3+1, 2):
                f = (j-3*d)//2
                w = ((d+1)//2, d//2, (f+2)//3, (f+1)//3, f//3)
                if wb[w] or wc[w]:
                    expected.add((y, d, w, wc[w], wb[w]))
        actual = {(r['length'], r['dotted'], tuple(r['weight']), r['cartan_rank'], r['orthogonal_upper'])
                  for r in records if r['level']==j}
        assert actual == expected
        verified.append(dict(level=j,orthogonal=sum(b.values()),symplectic=sum(c.values()),
                             all_nonempty_sector_lengths=[len(b),len(c)],
                             differing_sector_lengths=len(differences),component_pairs=sum(comps[0].values()),
                             multiplet_pairs=sum(highest(comps[0]).values()),central_weights=len(actual)))
    return verified


if __name__ == '__main__':
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    records = json.loads(gzip.decompress((HERE/'selected.json.gz').read_bytes()))
    report = dict(levels=check_census(records), primes=[1009, 1013], minors=[])
    for record in records:
        minor = matrix(record)
        report['minors'].append(minor)
        print(record['level'], record['length'], record['dotted'], minor['order'],minor['determinants'],flush=True)
    args.out.write_text(json.dumps(report, indent=2)+'\n')
