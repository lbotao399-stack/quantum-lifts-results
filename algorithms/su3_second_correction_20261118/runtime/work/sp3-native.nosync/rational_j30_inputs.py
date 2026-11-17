"""Exact rational inputs for the Sp(2) j=30 graviton boundary problem.

This extension does not mutate pinned projector, basis or evaluator sources.
Reconstructed vectors must pass exact raising equations and reduce to the
saved modular vectors after undoing the recorded integer denominator.
"""
import sp_environment
import sp_census as census
import hw_extremal as projection
from exact_trace_eval import rational_reconstruct
from collections import defaultdict, Counter
from fractions import Fraction
from functools import lru_cache
from math import lcm
from pathlib import Path
import hashlib
import json
import pickle
import time

HERE = Path(__file__).resolve().parent
OUT = HERE / 'sp2_j30_rational_inputs'
ROOT = HERE.parent / 'cluster-native.nosync/outgoing-package-check-h87hckd3/results'
PRIME = 1000003


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply(left, right):
    value = defaultdict(int)
    for (_, x), cx in left.items():
        for (_, y), cy in right.items():
            sign, word = census.base.canon((), x + y)
            if sign:
                value[word] += sign * cx * cy
    return {w: c for w, c in value.items() if c}


def reduce_expr(expr, pivots):
    value = {w: Fraction(c) for w, c in expr.items() if c}
    while value:
        word = min(value)
        if word not in pivots:
            break
        coefficient = value[word]
        for w, c in pivots[word].items():
            value[w] = value.get(w, 0) - coefficient * c
            if not value[w]:
                del value[w]
    return value


def main():
    start = time.monotonic()
    OUT.mkdir(exist_ok=True)
    family_root = HERE / 'sp2_5555_graviton_family'
    family = pickle.loads((family_root / 'family.pkl').read_bytes())
    rational = pickle.loads((family_root / 'rational_family.pkl').read_bytes())
    for name in ('family', 'rational_family'):
        metadata = json.loads((family_root / (name + '.json')).read_text())
        assert digest(family_root / (name + '.pkl')) == metadata[name + '_sha256']

    @lru_cache(None)
    def product(indices):
        if not indices:
            return {((), ()): 1}
        return multiply(product(indices[:-1]), family['generators'][indices[-1]]['expr'])

    expressions = []
    for row in rational['kernel']:
        value = defaultdict(Fraction)
        for monomial, c in zip(rational['products'], row):
            if c == '0':
                continue
            coefficient = Fraction(c)
            for word, factor in product(tuple(monomial)).items():
                value[word] += coefficient * factor
        value = {w: c for w, c in value.items() if c}
        assert not census.base.action_expr(value)
        assert not projection.action(value, 0, 1)
        assert not projection.action(value, 1, 2)
        expressions.append(value)
    pivots = {}
    independent = []
    for i, expr in enumerate(expressions):
        remainder = reduce_expr(expr, pivots)
        if remainder:
            word = min(remainder)
            coefficient = remainder[word]
            pivots[word] = {w: c / coefficient for w, c in remainder.items()}
            independent.append(i)
    (OUT / 'gravitons.pkl').write_bytes(pickle.dumps(dict(expressions=expressions,
        independent=independent), protocol=5))
    print(json.dumps(dict(stage='exact_gravitons', terms=list(map(len, expressions)),
        independent=independent, seconds=time.monotonic()-start)), flush=True)
    product.cache_clear()
    census.base.clear_caches()
    projection.word_action.cache_clear()

    source = pickle.loads((ROOT / 'basis6/basis.pkl').read_bytes())
    metadata = json.loads((ROOT / 'basis6/result.json').read_text())
    assert digest(ROOT / 'basis6/basis.pkl') == metadata['basis_sha256']
    residues = {c for expr in source for c in expr.values()}
    reconstructed = {c: rational_reconstruct(c, PRIME) for c in residues}
    basis = []
    denominators = []
    for i, expr in enumerate(source):
        denominator = lcm(*(reconstructed[c].denominator for c in expr.values()))
        value = {w: int(reconstructed[c] * denominator) for w, c in expr.items()}
        assert not projection.action(value, 0, 1)
        assert not projection.action(value, 1, 2)
        assert all((c - expr[w] * denominator) % PRIME == 0 for w, c in value.items())
        basis.append(value)
        denominators.append(denominator)
        if i % 100 == 0:
            print(json.dumps(dict(stage='rational_basis', rows=i+1,
                seconds=time.monotonic()-start)), flush=True)
    (OUT / 'basis6.pkl').write_bytes(pickle.dumps(dict(basis=basis,
        denominators=denominators, modular_basis_sha256=digest(ROOT/'basis6/basis.pkl')), protocol=5))
    record = dict(status='passed', field='Q', sector=[5,5,5,5],
        graviton_formal_rank=len(independent), graviton_independent_rows=independent,
        graviton_terms=list(map(len,expressions)), graviton_raising_checks=34,
        graviton_closed_checks=17, incoming_rational_basis_rows=len(basis),
        incoming_raising_checks=2*len(basis), incoming_modular_coefficients_checked=sum(map(len,basis)),
        incoming_degrees=dict(Counter(sum(census.base.count(next(iter(e)))[3:6]) for e in basis)),
        maximum_denominator=max(denominators),
        files_sha256={name:digest(OUT/name) for name in ('gravitons.pkl','basis6.pkl')},
        source_sha256=digest(Path(__file__)), seconds=time.monotonic()-start,
        scope='Exact rational input expressions; no rational boundary/product quotient yet.')
    (OUT/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)


if __name__ == '__main__':
    main()
