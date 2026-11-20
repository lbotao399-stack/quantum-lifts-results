"""An integral functional on the first equation that kills parents with two bilinear traces."""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from pathlib import Path

import content_parents
import even_content_frames as frames
import finite_trace_identity as ft
import parent_frontier

s = ft.src
TARGET = bytes.fromhex('0000000000000505050814')
U = (-1,1,0, -1,0,0, 0,0,1)
V = (0,0,0, 0,0,1, 0,0,0)
W = (0,0,0, 0,0,0, 1,0,0)
BASIS = (U, V, W)
MATRICES = {(0,-1): U, (8,-1): V, (20,-1): W,
            (5,0): U, (5,1): V, (5,2): W}


def multiply(a, b):
    return tuple(a[i]*b[j] + a[i+1]*b[3+j] + a[i+2]*b[6+j]
                 for i in (0,3,6) for j in (0,1,2))


@lru_cache(None)
def trace(tr, indices):
    it = iter(indices)
    result = (1,0,0, 0,1,0, 0,0,1)
    for i in tr:
        result = multiply(result, MATRICES[i, next(it) if ft.ODD[i] else -1])
    return result[0]+result[4]+result[8]


@lru_cache(None)
def evaluation(w):
    if bytes(sorted(b''.join(w))) != TARGET:
        return 0
    sizes = [tr.count(5) for tr in w]
    assert sum(sizes) == 3
    value = 0
    for p in permutations(range(3)):
        term = ft.sign(p)
        start = 0
        for tr, n in zip(w, sizes):
            term *= trace(tr, p[start:start+n])
            start += n
        value += term
    return value


def image(w):
    return sum(c * evaluation(z) for z,c in content_parents.f.q(w).items())


def run(rhs):
    gram = [[sum(a[3*i+j]*b[3*j+i] for i in range(3) for j in range(3))
             for b in BASIS] for a in BASIS]
    assert gram == [[0]*3 for _ in BASIS]
    expected = ((0,0,(0,)),(0,0,(1,2)),(0,1,(0,)),(1,0,(1,)))
    assert tuple(s.ALPH[i] for i in sorted(set(TARGET))) == expected
    m = {k: tuple(tuple(a[i:i+3]) for i in (0,3,6)) for k,a in MATRICES.items()}
    rhs_value = Fraction()
    rhs_words = Counter()
    h = hashlib.sha256()
    for w, c in s.stream_rows(rhs, h):
        if bytes(sorted(b''.join(w))) == TARGET:
            value = evaluation(w)
            assert value == ft.exterior_value(w, m, 3)
            rhs_value += Fraction(c) * value
            rhs_words[sum(len(t)==2 for t in w)] += 1
    assert h.hexdigest() == s.INPUT and rhs_value == 771120
    _, incidence = parent_frontier.incidence()
    counts = defaultdict(Counter)
    examples = {}
    enumerated = {}
    for g in sorted(incidence[TARGET]):
        words = content_parents.all_words(g)
        assert len(words) == frames.base.control(g)
        enumerated[g.hex()] = len(words)
        for w in words:
            sign, partner = frames.conjugate(w)
            if w > partner or (w == partner and sign < 0):
                continue
            value = image(w) + sign * image(partner)
            shape = tuple(sorted(map(len,w)))
            assert shape.count(2)<2 or value == 0
            counts[shape]['tested'] += 1
            if value:
                counts[shape]['nonzero'] += 1
                examples.setdefault(shape,[s.hexes(w),value])
    assert examples
    return dict(rhs_sha256=h.hexdigest(), target=TARGET.hex(), exterior_degree=3,
                colour_basis=BASIS, trace_gram=gram, assignments=[[i,j,a] for (i,j),a in sorted(MATRICES.items())],
                rhs_words_by_bilinear_count=dict(rhs_words), pairing_with_negative_Q1=str(rhs_value),
                unprojected_words_by_incident_content=enumerated,
                parity_even_orbits_by_trace_lengths=[dict(lengths=sh, tested=rec['tested'],
                    nonzero=rec['nonzero'], example=examples.get(sh)) for sh,rec in sorted(counts.items())],
                parent_projection='1+C')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz')
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();assert __debug__
    z=run(a.rhs)
    a.out.write_text(json.dumps(z,indent=2)+'\n')
    print(z['pairing_with_negative_Q1'],[(r['lengths'],r['nonzero']) for r in z['parity_even_orbits_by_trace_lengths']],flush=True)
