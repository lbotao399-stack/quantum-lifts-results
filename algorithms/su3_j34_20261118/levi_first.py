"""A shared SU(3) correction on the upper two-colour block plus a central scalar."""
import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, lcm
from pathlib import Path

import collision_solve as solve
from collision_first import block


@lru_cache(32000)
def trace_options(trace, scalar):
    positions = [i for i, k in enumerate(trace) if k == scalar]
    terms = Counter()
    for n in range(len(positions) + 1):
        for chosen in combinations(positions, n):
            rest = tuple(k for i, k in enumerate(trace) if i not in chosen)
            if rest:
                c, w = block.canonical_trace(rest)
            else:
                c, w = 2 + (-2)**n, ()
            if c:
                terms[n, w] += c
    return [(n, w, c) for (n, w), c in terms.items() if c]


def expansion(word, scalar):
    parts = [(0, (), 1)]
    for trace in word:
        parts = [(n + m, factors + ((z,) if z else ()), c * v)
                 for n, factors, c in parts for m, z, v in trace_options(trace, scalar)]
    for power, factors, c in parts:
        sign, w = block.canonical_multi(factors)
        if sign:
            yield power, w, c * sign


def extract(rhs, folder):
    folder.mkdir(exist_ok=True, parents=True)
    digest = hashlib.sha256()
    numerators = defaultdict(Counter)
    count = 0
    for line, numerator, w in block.encoded_stream(rhs):
        digest.update(line)
        count += 1
        for scalar in (1, 2, 4):
            for power, z, multiplicity in expansion(w, scalar):
                if power or scalar == 1:
                    numerators[scalar if power else 0, power][z] += numerator * multiplicity
    assert digest.hexdigest() == block.RHS_SHA and count == 517468
    cases = []
    for (scalar, power), row in sorted(numerators.items()):
        remaining = {w: c for w, c in row.items() if c}
        assert power % 2 == 0 or not remaining
        assert all(sum(map(len, w)) + power == 11 for w in remaining)
        name = f'levi-{scalar}-{power}.jsonl.gz'
        digest = block.dump_rows(folder / name, remaining)
        cases.append(dict(scalar=scalar, power=power, terms=len(remaining), file=name, sha256=digest))
    report = dict(input=block.RHS_SHA, input_words=count, denominator=120, cases=cases)
    (folder / 'extraction.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def domain(scalar, power):
    from g2_full_classical_screen_v1 import contents
    words = []
    for raw in contents((1, 1, 7, 4, 3), 10):
        t = sorted(map(block.code, raw))
        if t.count(scalar) < power:
            continue
        for _ in range(power):
            t.remove(scalar)
        words.extend(solve.wall.pair_words(tuple(t)))
    assert len(words) == len(set(words))
    return sorted(words)


def input_rows(folder, scalar, power):
    doc = json.loads((folder / 'extraction.json').read_bytes())
    assert doc['input'] == block.RHS_SHA and doc['denominator'] == 120
    raw = defaultdict(dict)
    for case in doc['cases']:
        if (scalar, power) != (case['scalar'], case['power']):
            continue
        stream = gzip.decompress((folder / case['file']).read_bytes())
        assert hashlib.sha256(stream).hexdigest() == case['sha256']
        for c, w in map(json.loads, stream.splitlines()):
            w = tuple(map(tuple, w))
            raw[solve.content(w)][w] = c
    return raw


def lift_pair_word(word, scalar, power):
    pairs = power // 2
    trace = (scalar, scalar)
    sign, full = block.canonical_multi(word + (trace,) * pairs)
    assert sign == 1
    return full, 6**pairs * comb(word.count(trace) + pairs, pairs)


def block_inputs(folder, scalar, power):
    raw = input_rows(folder, scalar, power)
    parents = domain(scalar, power)
    images = defaultdict(lambda: defaultdict(dict))
    for i, w in enumerate(parents):
        for z, c in solve.differential(w, -1).items():
            images[solve.content(z)][i][z] = c
    return raw, parents, images


def build(folder, scalar, power, out):
    raw, parents, images = block_inputs(folder, scalar, power)
    full_domain = {w: i for i, w in enumerate(domain(0, 0))}
    lifts = [lift_pair_word(w, scalar, power) for w in parents]
    assert all(w in full_domain for w, _ in lifts)
    eq = []
    for t in sorted(set(raw) | set(images)):
        diagonal, rows, rhs = solve.block_data(raw[t], images[t], t)
        keys = solve.compress(rows, rhs)
        for key in keys:
            terms = [[full_domain[lifts[i][0]], row[key] * lifts[i][1]]
                     for i, row in rows.items() if key in row]
            eq.append([t, diagonal, [str(k) for k in key], terms, rhs.get(key, 0)])
    doc = dict(scalar=scalar, power=power, upper_pair_words=len(parents), equations=eq)
    raw_bytes = (json.dumps(doc, separators=(',', ':')) + '\n').encode()
    out.write_bytes(gzip.compress(raw_bytes, mtime=0))
    print('upper plus central', scalar, power, len(parents), len(eq), flush=True)
    return len(eq)


def coupled(blocks, out):
    all_words = domain(0, 0)
    equations, labels, sizes = [], [], []
    for path in sorted(blocks):
        data = json.loads(gzip.decompress(path.read_bytes()))
        scalar, power = data['scalar'], data['power']
        sizes.append([scalar, power, len(data['equations'])])
        for t, diagonal, key, row, b in data['equations']:
            equations.append((dict(row), b))
            labels.append((scalar, power, t, diagonal, key))
    mode, coefficients, value = solve.eliminate(equations)
    result = dict(input=block.RHS_SHA, full_pair_words=len(all_words), selected_equations=sizes,
                  kind=mode)
    if mode == 'solution':
        result.update(selected_rank=value, primitive=[[all_words[i], str(c / 120)]
                                                     for i, c in sorted(coefficients.items())])
    else:
        result.update(pairing=str(value / 120), functional=[[labels[i], str(c)]
                                                             for i, c in sorted(coefficients.items())])
    out.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('coupled', mode, len(coefficients), value, flush=True)
    return result


def verify_case(folder, scalar, power, result, out):
    """Evaluate both sides at every polynomial coordinate on this scalar slice."""
    doc = json.loads(result.read_bytes())
    assert doc['input'] == block.RHS_SHA
    raw, parents, images = block_inputs(folder, scalar, power)
    blocks = sorted(set(raw) | set(images))
    receipt = dict(scalar=scalar, power=power, blocks=len(blocks))
    if doc['kind'] == 'solution':
        values = {tuple(map(tuple, w)): Fraction(c) * 120 for w, c in doc['primitive']}
        co = [factor * values.get(w, 0) for w, factor in
              (lift_pair_word(z, scalar, power) for z in parents)]
        den = lcm(*(Fraction(v).denominator for v in co))
        count = 0
        for t in blocks:
            r = Counter({w: -den * c for w, c in raw[t].items()})
            for i, small in images[t].items():
                if co[i]:
                    solve.add(r, small, int(co[i] * den))
            count += len(solve.sum_poly({w: c for w, c in r.items() if c}, solve.diagonal_for(t)))
        assert count == 0
        receipt['nonzero_residual_monomials'] = count
    else:
        full_domain = {w: i for i, w in enumerate(domain(0, 0))}
        covector = defaultdict(dict)
        for (s, p, t, diagonal, key), coefficient in doc['functional']:
            if (s, p) == (scalar, power):
                assert diagonal == solve.diagonal_for(tuple(t))
                covector[tuple(t)][tuple(map(int, key))] = Fraction(coefficient)
        values = Counter()
        rhs_value = Fraction()
        for t, functional in covector.items():
            _, rows, rhs = solve.block_data(raw[t], images[t], t)
            rhs_value += sum(rhs.get(k, 0)*c for k, c in functional.items()) / 120
            for i, row in rows.items():
                w, factor = lift_pair_word(parents[i], scalar, power)
                values[full_domain[w]] += factor * sum(row.get(k, 0)*c for k, c in functional.items())
        receipt.update(parent_pairings=[[i, str(c)] for i, c in sorted(values.items()) if c],
                       rhs_pairing=str(rhs_value), covector_contents=len(covector))
    out.write_text(json.dumps(receipt, separators=(',', ':')) + '\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    p = sub.add_parser('extract')
    p.add_argument('--rhs', type=Path, required=True)
    p.add_argument('--folder', type=Path, required=True)
    p = sub.add_parser('coupled')
    p.add_argument('--blocks', type=Path, nargs='+', required=True)
    p.add_argument('--out', type=Path, required=True)
    for mode in ('build', 'verify'):
        p = sub.add_parser(mode)
        p.add_argument('--folder', type=Path, required=True)
        p.add_argument('--scalar', type=int, required=True)
        p.add_argument('--power', type=int, required=True)
        p.add_argument('--out', type=Path, required=True)
        if mode == 'verify':
            p.add_argument('--result', type=Path, required=True)
    assert __debug__
    args = vars(parser.parse_args())
    {'extract': extract, 'build': build, 'coupled': coupled, 'verify': verify_case}[args.pop('mode')](**args)
