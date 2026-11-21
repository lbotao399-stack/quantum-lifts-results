"""Restrict the parent first equation to a repeated scalar eigenvalue.

The chosen scalar equals x diag(1,1,-2). Every other Taylor coefficient
is an independent traceless two-colour matrix in the upper block.
The coefficients of distinct powers of the central scalar stay separate.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

import block_first as block
import block_reduce as reduce


@lru_cache(None)
def restrict_trace(trace, scalar):
    rest = tuple(k for k in trace if k != scalar)
    power = len(trace) - len(rest)
    if not rest:
        return power, 2 + (-2) ** power, ()
    sign, reduced = block.canonical_trace(rest)
    return power, sign, reduced


def restrict_word(word, scalar):
    power, coefficient, traces = 0, 1, []
    for trace in word:
        n, sign, reduced = restrict_trace(trace, scalar)
        power += n
        coefficient *= sign
        if reduced:
            traces.append(reduced)
    sign, reduced = block.canonical_multi(tuple(traces))
    return power, sign * coefficient, reduced


def extract(rhs, out, *, rhs_sha=block.RHS_SHA, denominator=120, sign_of_input=1,
            input_words=517468):
    if sign_of_input not in (-1, 1):
        raise ValueError("the input sign must be +1 or -1")
    out.mkdir(exist_ok=True, parents=True)
    digest = hashlib.sha256()
    scalars = [block.code((0, 0, (i,))) for i in range(3)]
    numerators = defaultdict(Counter)
    input_counts = Counter()
    for line, coefficient, word in block.encoded_stream(rhs, denominator):
        digest.update(line)
        for scalar in scalars:
            power, sign, reduced = restrict_word(word, scalar)
            key = (scalar, power)
            input_counts[key] += 1
            if sign:
                numerators[key][reduced] += coefficient * sign_of_input * sign
    assert digest.hexdigest() == rhs_sha
    records = []
    for (scalar, power), count in sorted(input_counts.items()):
        row = {w: c for w, c in numerators[scalar, power].items() if c}
        assert all(sum(map(len, w)) == 11 - power for w in row)
        assert power % 2 == 0 or not row
        name = f'scalar{scalar}-power{power}.jsonl.gz'
        sha = block.dump_rows(out / name, row)
        records.append(dict(scalar=scalar, power=power, input_words=count,
                            restricted_trace_words=len(row), file=name, trace_sha256=sha))
    assert all(sum(r['input_words'] for r in records if r['scalar'] == k) == input_words for k in scalars)
    report = dict(rhs_sha256=digest.hexdigest(), denominator=denominator, rows=records)
    if sign_of_input != 1:
        report["sign_applied_to_input"] = sign_of_input
    (out / 'extraction.json').write_text(json.dumps(report, indent=2) + '\n')
    print([(r['scalar'],r['power'],r['restricted_trace_words']) for r in records], flush=True)
    return report


@lru_cache(None)
def pair_words(content):
    """Polarized SO(3) invariant spanning set in even degree."""
    if not content:
        return ((),)
    first, remaining = content[0], content[1:]
    words = set()
    for second in sorted(set(remaining)):
        others = list(remaining)
        others.remove(second)
        sign, pair = block.canonical_trace((first, second))
        if sign:
            for tail in pair_words(tuple(others)):
                sign, word = block.canonical_multi((pair,) + tail)
                if sign:
                    words.add(word)
    return tuple(sorted(words))


def parent_words(scalar, power, mu=(1, 1, 7, 4, 3)):
    from g2_full_classical_screen_v1 import contents
    assert power % 2 == 0
    words = []
    parts = []
    for raw in contents(mu, 10):
        content = sorted(map(block.code, raw))
        if content.count(scalar) != power:
            continue
        reduced = tuple(k for k in content if k != scalar)
        part = pair_words(reduced)
        words.extend(part)
        parts.append((reduced, len(part)))
    assert len(words) == len(set(words))
    return sorted(words), sorted(parts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rhs', type=Path, default=block.HERE / 'negative_Q1.jsonl.gz')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert __debug__
    extract(args.rhs, args.out)
