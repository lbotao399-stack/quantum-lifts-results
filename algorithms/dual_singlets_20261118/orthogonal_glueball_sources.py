"""Finite-colour trace relations in the possible sources of a pure glueball power."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import comb
from pathlib import Path

from b2_schur_complex import centralizer, character, multiplicity, partitions, transpose


@lru_cache(None)
def missing_diagrams(degree, colours):
    # Invariants of 2*degree vector slots have even rows. Transposition
    # shortens the diagrams that violate the defining dimension bound.
    return tuple(transpose(tuple(2*k for k in reversed(parts)))
                 for parts in partitions(degree) if len(parts) > colours)


@lru_cache(None)
def missing_vector_trace(degree, colours, cycles):
    sign = (-1)**(2*degree-len(cycles))
    return sign * sum(character(shape, cycles) for shape in missing_diagrams(degree, colours))


@lru_cache(None)
def missing_adjoint_trace(cycles, colours):
    degree = sum(cycles)
    if not missing_diagrams(degree, colours):
        return 0
    types = sorted(Counter(cycles).items())
    total = 0
    for crossed in product(*(range(count+1) for _,count in types)):
        vector_cycles = []
        coefficient = (-1)**sum(crossed)
        for (length, count), n in zip(types, crossed):
            vector_cycles.extend([length]*(2*(count-n)) + [2*length]*n)
            coefficient *= comb(count, n)
        total += coefficient * missing_vector_trace(degree, colours, tuple(sorted(vector_cycles)))
    value, remainder = divmod(total, 2**len(cycles))
    if remainder:
        raise ArithmeticError('Nonintegral character of the adjoint relation space')
    return value


@lru_cache(None)
def relations(shape, length, colours):
    return multiplicity(shape, length, lambda cycles: missing_adjoint_trace(cycles, colours))


def physical_image_relations(power, order, colours):
    """The full spin/flavour singlet of first jets after applying the tree charge."""
    ghosts = (2,)*(power-order)
    scalars = (order,)*3
    total = Fraction()
    for gamma in partitions(sum(ghosts)):
        first = character(ghosts, gamma)
        if not first:
            continue
        for phi in partitions(sum(scalars)):
            second = character(scalars, phi)
            if second:
                trace = missing_adjoint_trace(tuple(sorted(gamma+phi)), colours)
                total += Fraction(first*second*trace, centralizer(gamma)*centralizer(phi))
    if total.denominator != 1 or total < 0:
        raise ArithmeticError('Invalid physical relation multiplicity')
    return int(total)


def proper_parents(shape, colours):
    """All proper Taylor subdiagrams in the finite-identity range."""
    last = sum(shape)
    rows = []
    for degree in range(colours+1, last):
        for parts in partitions(degree):
            parent = tuple(reversed(parts))
            if len(parent) > len(shape) or any(i>j for i,j in zip(parent,shape)):
                continue
            count = relations(parent, degree, colours)
            rows.append(dict(shape=parent, degree=degree, image_space_relations=count))
    return rows


def blocks(colours, power, through):
    rows = []
    for order in range(1, min(power-3, through)+1):
        diagram = (power-order, power-order) + (3,)*order
        total = sum(diagram)
        assert total == 2*power+order
        target = relations(diagram, total, colours)
        physical = physical_image_relations(power, order, colours)
        assert physical >= target
        source = relations(diagram, total-1, colours)
        row = dict(order=order, diagram=diagram, length=total-1, source_relations=source,
                   image_space_relations=target, full_singlet_image_space_relations=physical)
        if total % 2 and target > source and colours == power+2*order+2:
            row['fortuitous_multiplicity_bounds'] = [target-source, target]
            row['proper_parents'] = proper_parents(diagram, colours)
            assert all(p['image_space_relations']==0 for p in row['proper_parents'])
        rows.append(row)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', nargs=2, type=int, action='append', required=True,
                        metavar=('COLOURS','POWER'))
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--through', type=int, default=3)
    args = parser.parse_args()
    if not __debug__ or args.through < 1 or any(n < 3 or k < 4 for n,k in args.case):
        parser.error('Use at least three colours and power four, with checks enabled')
    results = {}
    for colours, power in args.case:
        result = blocks(colours, power, args.through)
        print(colours, power, [(r['source_relations'],r['image_space_relations'],
                               r['full_singlet_image_space_relations']) for r in result], flush=True)
        results[f'O{colours}_S{power}'] = result
    args.out.write_text(json.dumps(results, indent=2)+'\n')
