#!/usr/bin/env python3
"""Exact additional-sector scan of the two/three-trace coefficient pattern.

This is a coefficient identity experiment, not a proof of an all-rank theorem.
It uses the previously fixed auxiliary and Q0 conventions unchanged.
"""
import json
from fractions import Fraction
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'scalar'))
import minimal_scalar_certificates as minimal
import pfaffian_scalar_candidates as p

started = time.monotonic()
records = []
for shape in [(3,3,3),(5,3,3),(7,3,3),(9,3,3),(11,3,3),(5,5,3),(7,5,3),(9,5,3),(7,7,3),(5,5,5),(7,5,5),(9,5,5),(7,7,5),(7,7,7)]:
    sources, terminal, _, _ = minimal.operators(shape)
    sources = list(map(p.Coefficient, sources))
    terminal = p.Coefficient(terminal)
    rows = []
    for triples in range(1, min(shape) + 1, 2):
        word = [(0,1,2)] * triples
        word += [(i,i) for i,n in enumerate(shape) for _ in range((n-triples)//2)]
        sign, word = p.canonical(word)
        assert sign == 1
        target = terminal(word)
        image = [p.image(source,word) for source in sources]
        assert target
        rows.append({'ABC_trace_factors': triples, 'word': p.encode(word), 'target': str(target), 'Q0_sources': list(map(str,image)), 'Q0_divided_by_ABC_count_and_target': [str(v/(triples*target)) for v in image]})
    assert len({tuple(row['Q0_divided_by_ABC_count_and_target']) for row in rows}) == 1
    record = {'shape': shape, 'orthogonal_colours': sum(shape)-2, 'rows': rows, 'observed_identity_exact': True}
    records.append(record)
    (HERE/'PROBE_LAW_SCAN.json').write_text(json.dumps({'status':'finite exact experiment; no all-rank claim','records':records,'seconds':time.monotonic()-started},indent=2)+'\n')
    print(shape, 'K=', str(3*Fraction(rows[1]['target'])/Fraction(rows[0]['target'])), 'rows=',len(rows), 'seconds=', time.monotonic()-started, flush=True)
