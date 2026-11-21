"""Independently multiply the first corrections on generic upper traceless matrices."""
import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import lcm
from pathlib import Path

from fortuity.full.letters import Q_letter
import block_first as block
import block_reduce as physical


def derivative(w, scalar):
    result = Counter()
    sign = 1
    for i, trace in enumerate(w):
        for j, k in enumerate(trace):
            central = Counter()
            for c, a, b in Q_letter(block.letter(k)):
                pair = block.code(a), block.code(b)
                if scalar in pair:
                    central[tuple(v for v in pair if v != scalar)] += c
                else:
                    changed = trace[:j] + pair + trace[j+1:]
                    result[w[:i] + (changed,) + w[i+1:]] += sign * c
            assert not any(central.values())
            if block.odd(k):
                sign = -sign
    return result


def run(extraction, files, *, rhs_sha=block.RHS_SHA, denominator=120):
    doc = json.loads((extraction/'extraction.json').read_bytes())
    assert doc['rhs_sha256'] == rhs_sha and doc['denominator'] == denominator
    results = []
    for case in doc['rows']:
        scalar, power = case['scalar'], case['power']
        if not case['restricted_trace_words']:
            continue
        record = json.loads((files/f'case-{scalar}-{power}.json').read_bytes())
        assert record['input'] == case and record['kind'] == 'solution'
        b = gzip.decompress((extraction/case['file']).read_bytes())
        assert hashlib.sha256(b).hexdigest() == case['trace_sha256']
        source = {tuple(map(tuple,w)): c for c,w in map(json.loads,b.splitlines())}
        sol = [(tuple(map(tuple,w)), Fraction(c)) for w,c in record['primitive']]
        den = lcm(denominator, *(c.denominator for w,c in sol))
        rows = defaultdict(Counter)
        for w,c in source.items():
            rows[tuple(sorted(sum(w, ())))][w] -= c * (den//denominator)
        for w,c in sol:
            for z,a in derivative(w, scalar).items():
                rows[tuple(sorted(sum(z, ())))][z] += a * int(c*den)
        words = 0
        for row in rows.values():
            remainder = {w:c for w,c in row.items() if c}
            words += len(remainder)
            assert not physical.polynomial(remainder)
        results.append(dict(scalar=scalar, power=power, primitive_terms=len(sol),
                            raw_trace_terms_tested=words, residual_monomials=0))
    return dict(rhs=rhs_sha, generic_traceless_upper_block=True, cases=results)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('extraction','files','out'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();assert __debug__
    z=run(a.extraction,a.files)
    a.out.write_text(json.dumps(z,indent=2)+'\n')
    print(z['cases'],flush=True)
