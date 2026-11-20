"""Coupled rational first equation on a two-colour block at a scalar collision."""
import argparse
import gzip
import hashlib
import heapq
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, lcm
from pathlib import Path

import numpy as np
from flint import nmod_mat
from fortuity.full.letters import Q_letter
from fortuity.full.words import Q_multi

import collision_first as wall
from block_first import canonical_multi, canonical_trace, code, letter, odd
from verify import exterior_sign, multiply


def content(word):
    return tuple(sorted(k for trace in word for k in trace))


def differential(word, scalar):
    """Project the full signed classical variation onto the upper traceless block."""
    result = Counter()
    for full, c in Q_multi(wall.reduce.to_fort(word)).items():
        raw = wall.reduce.from_fort(full)
        if any(scalar in trace for trace in raw):
            continue
        sign, reduced = wall.reduce.canonical(raw)
        if sign:
            result[reduced] += c * sign
    return {w: c for w, c in result.items() if c}


def split_controls(scalar):
    """A split involving the central scalar is a commuting matrix commutator."""
    from fortuity.full.letters import letters_up_to
    terms = 0
    for L in letters_up_to(1, 1):
        pairs = {(a, b): c for c, a, b in Q_letter(L)}
        for (a, b), c in pairs.items():
            assert pairs[b, a] == -c * (-1) ** (odd(code(a)) * odd(code(b)))
            terms += int(code(a) == scalar or code(b) == scalar)
    return terms


@lru_cache(20000)
def trace_polynomial(trace, diagonal):
    total = Counter()
    for path in product((0, 1), repeat=len(trace)):
        even = mask = 0
        sign = 1
        for k, i, j in zip(trace, path, path[1:] + path[:1]):
            if k == diagonal and i != j:
                break
            if i == j == 1:
                sign = -sign
            if k == diagonal:
                continue
            slot = 3 * k + (0 if i == j else 1 + i)
            if odd(k):
                bit = 1 << slot
                sign *= exterior_sign(mask, bit)
                if not sign:
                    break
                mask |= bit
            else:
                even += 1 << (4 * slot)
        else:
            total[even, mask] += sign
    return {k: v for k, v in total.items() if v}


def polynomial(word, diagonal):
    value = {(0, 0): 1}
    for trace in word:
        value = multiply(value, trace_polynomial(trace, diagonal))
        if not value:
            break
    return value


def add(out, row, factor):
    for k, v in row.items():
        out[k] += v * factor


def sum_poly(words, diagonal):
    out = Counter()
    for w, c in words.items():
        add(out, polynomial(w, diagonal), c)
    return {k: v for k, v in out.items() if v}


def diagonal_for(t):
    counts = Counter(t)
    return min((k for k in t if not odd(k)), key=lambda k: (-counts[k], k), default=-1)


def compress(rows, rhs):
    """Select monomial equations. An accepted solution must pass every omitted equation."""
    keys = sorted(set(rhs).union(*(row.keys() for row in rows.values())))
    if not keys:
        return []
    prime = 1000003
    key_id = {k: i for i, k in enumerate(keys)}
    matrix = np.zeros((len(rows) + 1, len(keys)), dtype=np.int64)
    for i, poly in enumerate([*rows.values(), rhs]):
        for k, v in poly.items():
            matrix[i, key_id[k]] = v % prime
    echelon, rank = nmod_mat(matrix.tolist(), prime).rref()
    pivots = []
    j = 0
    for i in range(rank):
        while not echelon[i, j]:
            j += 1
        pivots.append(keys[j])
        j += 1
    return pivots


def block_data(raw, images, t):
    diagonal = diagonal_for(t)
    rhs = sum_poly(raw, diagonal)
    small = {w: polynomial(w, diagonal) for row in images.values() for w in row}
    rows = {}
    for i, row in images.items():
        value = Counter()
        for w, c in row.items():
            add(value, small[w], c)
        nz = {k: v for k, v in value.items() if v}
        if nz:
            rows[i] = nz
    return diagonal, rows, rhs


def primitive_row(row, rhs):
    scale = abs(rhs)
    for x in row.values():
        scale = gcd(scale, x)
        if scale == 1:
            break
    scale = scale or 1
    return {i: v // scale for i, v in row.items() if v}, rhs // scale, scale


def eliminate(equations):
    """Integer elimination with a sparsest-variable pivot and reversible row operations."""
    rows = {}
    support = defaultdict(set)
    operations = []
    back = []
    contradiction = None
    for i, (row, b) in enumerate(equations):
        rows[i] = (dict(row), b)
        for k in row:
            support[k].add(i)
        if not row and b:
            contradiction = i
            break
    queue = [(len(ids), k) for k, ids in support.items() if ids]
    heapq.heapify(queue)
    while queue and contradiction is None:
        n, k = heapq.heappop(queue)
        if n != len(support[k]):
            continue
        pivot = min(support[k], key=lambda i: (len(rows[i][0]), abs(rows[i][0][k]), i))
        p, y = rows.pop(pivot)
        back.append((k, p, y))
        affected = set(p)
        for j in p:
            support[j].remove(pivot)
        for i in sorted(support[k]):
            row, b = rows[i]
            common = gcd(p[k], row[k])
            u, v = p[k] // common, -row[k] // common
            updated = {j: u * x for j, x in row.items()}
            for j, x in p.items():
                updated[j] = updated.get(j, 0) + v * x
            new, b, scale = primitive_row(updated, u * b + v * y)
            operations.append((i, pivot, u, v, scale))
            removed, added = row.keys() - new.keys(), new.keys() - row.keys()
            for j in removed:
                support[j].remove(i)
            for j in added:
                support[j].add(i)
            affected.update(removed | added)
            rows[i] = new, b
            if not new and b:
                contradiction = i
                break
        for j in affected:
            if support[j]:
                heapq.heappush(queue, (len(support[j]), j))
    if contradiction is not None:
        dual = {contradiction: Fraction(1)}
        for i, p, u, v, scale in reversed(operations):
            if i not in dual:
                continue
            factor = dual[i] / scale
            dual[i] = factor * u
            dual[p] = dual.get(p, 0) + factor * v
            if not dual[p]:
                del dual[p]
        check = Counter()
        value = 0
        for i, factor in dual.items():
            row, b = equations[i]
            add(check, row, factor)
            value += b * factor
        assert not any(check.values()) and value
        return 'separator', dual, value
    solution = {}
    for k, row, b in reversed(back):
        value = Fraction(b - sum(c * solution.get(j, 0) for j, c in row.items()), row[k])
        if value:
            solution[k] = value
    assert all(sum(c * solution.get(j, 0) for j, c in row.items()) == b for row, b in equations)
    return 'solution', solution, len(back)


def solve(extraction, scalar, power, out):
    report = json.loads((extraction / 'extraction.json').read_bytes())
    rec, = [r for r in report['rows'] if (r['scalar'], r['power']) == (scalar, power)]
    assert report['rhs_sha256'] == wall.block.RHS_SHA and report['denominator'] == 120
    data = gzip.decompress((extraction / rec['file']).read_bytes())
    assert hashlib.sha256(data).hexdigest() == rec['trace_sha256']
    source = defaultdict(dict)
    for c, w in map(json.loads, data.splitlines()):
        word = tuple(map(tuple, w))
        source[content(word)][word] = c
    domain, parts = wall.parent_words(scalar, power)
    images = defaultdict(lambda: defaultdict(dict))
    for i, word in enumerate(domain):
        for z, c in differential(word, scalar).items():
            images[content(z)][i][z] = c
    targets = sorted(set(source) | set(images))
    equations, labels, statistics = [], [], []
    for t in targets:
        diagonal, rows, rhs = block_data(source[t], images[t], t)
        keys = compress(rows, rhs)
        for a, m in keys:
            key = a, m
            equations.append(({i: row[key] for i, row in rows.items() if key in row}, rhs.get(key, 0)))
            labels.append((t, diagonal, a, m))
        statistics.append([t, len(rows), len(keys), len(rhs)])
    print('matrix', scalar, power, len(domain), len(equations), flush=True)
    mode, co, value = eliminate(equations)
    doc = dict(scalar=scalar, power=power, input=rec, split_terms_with_scalar=split_controls(scalar),
               parent_contents=len(parts), parent_pair_words=len(domain), target_blocks=statistics,
               equations_used=len(equations), kind=mode)
    if mode == 'solution':
        denominator = lcm(*(v.denominator for v in co.values()))
        nonzero = checked = 0
        for t in targets:
            row = Counter({w: -denominator * c for w, c in source[t].items()})
            for i, small in images[t].items():
                if i in co:
                    add(row, small, int(denominator * co[i]))
            remainder = {w: c for w, c in row.items() if c}
            residual = sum_poly(remainder, diagonal_for(t))
            nonzero += len(residual)
            checked += 1
        assert nonzero == 0
        doc.update(rank_of_selected_equations=value, checked_blocks=checked, nonzero_residual_monomials=0,
                   primitive=[[domain[i], str(v / 120)] for i, v in sorted(co.items())])
    else:
        doc.update(pairing_with_negative_Q1=str(value / 120),
                   functional=[[labels[i], str(v)] for i, v in sorted(co.items())])
    out.write_text(json.dumps(doc, separators=(',', ':')) + '\n')
    print(scalar, power, mode, len(co), 'value', value, flush=True)
    return doc


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extraction', type=Path, required=True)
    p.add_argument('--scalar', type=int, required=True)
    p.add_argument('--power', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    assert __debug__
    solve(**vars(p.parse_args()))
