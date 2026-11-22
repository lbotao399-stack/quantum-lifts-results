"""Cartan restrictions of the first-jet intermediate targets and extremal relation counts."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

from b2_schur_complex import centralizer, character, partitions, transpose


def merge(left, right):
    counts = Counter()
    for a, c in left.items():
        for b, d in right.items():
            counts[tuple(sorted(a+b))] += c*d
    return counts


@lru_cache(None)
def identical_blocks(size, count):
    """Cycle index for permutations of count unordered blocks of equal size."""
    terms = Counter()
    for outer in partitions(count):
        part = {(): Fraction(1, centralizer(outer))}
        for k in outer:
            inner = {tuple(k*j for j in cycle):Fraction(1,centralizer(cycle))
                     for cycle in partitions(size)}
            part = merge(part, inner)
        terms.update(part)
    assert sum(terms.values()) == 1
    return terms


def terminal(power, order, rank):
    shape = (power-order,)*2 + (3,)*order
    colour = transpose(shape)
    contributions = []
    for parts in partitions(sum(shape)):
        if any(i%2 or i>colour[0] for i in parts):
            continue
        cycles = {(): Fraction(1)}
        for size, count in Counter(parts).items():
            cycles = merge(cycles, identical_blocks(size, count))
        dim = sum(v*character(colour, w) for w,v in cycles.items())
        assert dim >= 0 and dim.denominator == 1
        if dim:
            contributions.append(dict(lengths=parts, multiplicity=int(dim),
                                      independent_on_Cartan=len(parts)<=rank))
    all_dims = sum(r['multiplicity'] for r in contributions)
    cartan = sum(r['multiplicity'] for r in contributions if r['independent_on_Cartan'])
    return dict(grade=order, diagram=shape, stable_classical_dimension=all_dims,
                Cartan_rank=cartan, restriction_kernel=all_dims-cartan,
                trace_products=contributions)


def two_spins(degree, weight):
    # Character of the symmetric power of (3,3) plus a scalar for SL(2) x SL(2).
    coefficients = [Counter() for _ in range(degree+1)]
    coefficients[0][0,0] = 1
    for i,j in [*product((-2,0,2), repeat=2),(0,0)]:
        for n in range(1,degree+1):
            for (u,v),k in coefficients[n-1].items():
                coefficients[n][u+i,v+j] += k
    i,j = weight
    w = coefficients[degree]
    return w[i,j]-w[i+2,j]-w[i,j+2]+w[i+2,j+2]


def validate(inputs):
    results = {}
    for case, rows in inputs.items():
        colours, power = map(int,case[1:].split('_S'))
        checks = {}
        for r in rows:
            order = r['order']
            if colours != power+2*order+2:
                continue
            t = power-order-3
            expect = [two_spins(t-1,(2,2)) if t else 0, two_spins(t,(0,0))]
            assert expect == [r['source_relations'],r['image_space_relations']]
            checks[str(order)] = expect
        if checks:
            results[case] = checks
    return results


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('relations',type=Path)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    inputs=json.loads(args.relations.read_bytes())
    checks=validate(inputs)
    cases={}
    for colours, power in ((17,9),(18,10)):
        assert all(r['full_singlet_image_space_relations']==0 for r in inputs[f'O{colours}_S{power}'][:2])
        cases[f'O{colours}_S{power}']=[terminal(power,order,colours//2) for order in (0,1,2)]
    result=dict(extremal_colour_checks=checks, intermediate_targets=cases)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    for k, rows in cases.items():
        print(k,[(r['stable_classical_dimension'],r['Cartan_rank']) for r in rows])
