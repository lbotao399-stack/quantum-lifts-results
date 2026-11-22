"""Test the remaining SO(13) second pairing on two transported relation spaces."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import gzip
from itertools import combinations, product
import json

from verify import (A, ALPHABET, GHOSTS, HERE, determinant, encode, evaluate, first_duals,
                    insertion_parents, path_types, pullback, so_canon, unit, words)
from flint import fmpq_mat
from rank_four_lowered_replay import derivative
from terminal_commutators import parity


def reverse_odd(row, spin):
    """Transpose a native dotted-to-flavour-one derivation, retaining cyclic multiplicities."""
    result = Counter()
    for w,c in row.items():
        parents = set()
        for i,tr in enumerate(w):
            for j,x in enumerate(tr):
                if x == A:
                    sign,v = so_canon(w[:i]+(tr[:j]+(GHOSTS[spin],)+tr[j+1:],)+w[i+1:])
                    if sign:
                        parents.add(v)
        for v in parents:
            result[v] += c*derivative(v,spin,0).get(w,0)
    return {v:c for v,c in result.items() if c}


def pure_specs():
    for paths,p,m,_ in path_types(6,6,6,4):
        closed = [v for c,v in paths if c]
        opens = [v for c,v in paths if not c]
        assert closed in ([],['mp'])
        for axes in product(*(range(3) if len(v)%2 else [3] for v in opens)):
            for scalar in product(range(3),repeat=sum(v.count('s') for v in opens)):
                slots = iter(scalar)
                directions = [[unit(3,next(slots)) for _ in range(v.count('s'))] for v in opens]
                for plus,minus in product(combinations(range(3),p),combinations(range(3),m)):
                    yield dict(paths=opens,closed=closed,matrices=[unit(4,i) for i in axes],
                               single=[plus,minus],scalar_directions=directions)


def pure_relation(spec):
    rel = determinant(8-len(spec['closed']),spec)
    rel.closed = bool(spec['closed'])
    assert rel.pf_degree > 13
    return rel


def make_second():
    first,_ = first_duals(13,(7,7,2,0,0),[0,2])
    f = Counter({w:18*c for w,c in first[0].items()})
    f.update(first[1])
    f = {w:c for w,c in f.items() if c}
    row,info = pullback(13,f,-1)
    assert len(f) == 71 and len(row) == 2809 and info['checked_columns'] == 15928
    return row


def compact_poly(rel, basis):
    return {w:c for w in basis if (c := rel.commuting_coefficient(w))}


def substitution_on_section(poly, spin=None, letter=None):
    """Act on stable commuting traces with the same odd or even letter substitution."""
    result = Counter()
    for w,c in poly.items():
        prefix = 0
        for i,t in enumerate(w):
            target = list(t)
            if spin is None:
                k = (2,3,4,0,1)[letter]
                factor = target[k]
                target[2] += 2
            else:
                factor = target[spin]*(-1)**(prefix+sum(t[:spin]))
                target[spin] -= 1
                target[2] += 1
            prefix += parity(t)
            if not factor:
                continue
            sequence = w[:i]+(tuple(target),)+w[i+1:]
            key = lambda v:(sum(v),v)
            odd = [key(v) for v in sequence if parity(v)]
            if len(odd) != len(set(odd)):
                continue
            sign = (-1)**sum(a>b for j,a in enumerate(odd) for b in odd[j+1:])
            result[tuple(sorted(sequence,key=key))] += sign*factor*c
    return {w:c for w,c in result.items() if c}


def check_pure(task):
    spec, encoded = task
    rel = pure_relation(spec)
    score = sum(c*rel.coefficient(w) for w,c in encoded)
    sites = words((6,6,2,2,2))
    cartan = [rel.cartan(w) for w in sites[:3]]
    rel.clear()
    return [str(v) for v in cartan+[score]]


def test(workers, rebuild):
    row = make_second()
    twice = reverse_odd(reverse_odd(row,0),1)
    encoded = [(encode(w),c) for w,c in twice.items()]
    assert len(encoded) == 2574
    specs = list(pure_specs())
    saved = json.loads((HERE/'second.json').read_text())
    source_rows = saved['pure_eighth_rows']
    assert len(specs) == len(source_rows) == 56
    pivots = [0,2,28]
    checks = range(56) if rebuild else pivots
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,values in zip(checks,pool.map(check_pure,[(specs[i],encoded) for i in checks]),strict=True):
            assert values == source_rows[i]
    for c0,c1,c2,score in source_rows:
        # The pure second dual has the two displayed commuting coefficients on all 301 relations.
        pure = -Fraction(9,1280)*Fraction(c1)-Fraction(1,80)*Fraction(c2)
        assert Fraction(score) == -144*pure
    assert fmpq_mat([source_rows[i][:3] for i in pivots]).rank() == 3

    new_basis = words((5,5,4,2,2))
    sections, values = [], []
    pure7 = json.loads(gzip.decompress((HERE.parent/'so11_j36'/'relations13.json.gz').read_bytes()))
    forms = []
    for letter in ALPHABET:
        form = Counter()
        for w,c in row.items():
            for v,sg in insertion_parents(w,letter,'middle'):
                form[v] += sg*c
        forms.append({w:c for w,c in form.items() if c})
    for i,record in enumerate(pure7['independent']):
        rel = determinant(7,record['spec'])
        poly = compact_poly(rel,words((5,5,2,2,2)))
        for letter,form in enumerate(forms):
            value = evaluate(form,lambda w:rel.coefficient(encode(w)))
            assert str(value) == saved['pure_seventh_middle'][i][letter]
            # The end insertion is twice the middle insertion on commuting traces.
            # Verify its score independently for every non-scalar-A letter.
            if letter:
                end = sum(c*sum(s*rel.coefficient(encode(v)) for v,s in
                                insertion_parents(w,ALPHABET[letter],'ends')) for w,c in row.items())
                assert end == 2*value
            section = substitution_on_section(poly,letter=letter)
            sections.append([str(section.get(w,0)) for w in new_basis])
            values.append(str(value))
        rel.clear()

    image_sections = []
    for i in pivots:
        rel = pure_relation(specs[i])
        poly = compact_poly(rel,words((6,6,2,2,2)))
        section = substitution_on_section(substitution_on_section(poly,spin=1),spin=0)
        coord = [str(section.get(w,0)) for w in new_basis]
        image_sections.append(coord)
        sections.append(coord)
        values.append(source_rows[i][-1])
        rel.clear()
    minor = fmpq_mat([[r[j] for j in (0,2,10)] for r in image_sections])
    assert minor.det() != 0
    rank = fmpq_mat(sections).rank()
    augmented = fmpq_mat([r+[c] for r,c in zip(sections,values)]).rank()
    assert rank == augmented == 8
    print('mixed dual 15928 columns, pure polarizations 56, re-evaluated',len(checks),
          'joint ranks',rank,augmented,'image minor',minor.det())
    return sections,values,row


if __name__ == '__main__':
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('--rebuild',action='store_true',help='regenerate all 56 transported pairings')
    args.add_argument('--workers',type=int,default=2)
    options = args.parse_args()
    if not __debug__ or options.workers < 1:
        args.error('Run with assertions and a positive worker count')
    test(options.workers,options.rebuild)
