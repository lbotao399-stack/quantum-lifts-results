"""Specialize the complete reduced scalar-commutator space in two orthogonal Levis."""
import argparse
from collections import Counter
from itertools import product
import json
from math import gcd
from pathlib import Path
import sys

from flint import fmpq_mat, fmpz_mat
import numpy as np

from trace_ideal import abelian, order, trace

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import mixed_incoming as native
from f4_subgroup_equation import dimension, integral_rows
from orthogonal_glueball_sources import physical_image_relations


TO_LABEL = {(1, 0, ()): 3, (0, 1, ()): 4, **{(0, 0, (i,)): i for i in range(3)}}


def reduce_words(terms):
    result = Counter()
    for w, v in terms:
        signs_and_words = [trace(f) for f in w]
        s, m = order(tuple(f for _, f in signs_and_words))
        for phase, _ in signs_and_words:
            s *= phase
        if s:
            result[m] += s * v
    return {k: v for k, v in result.items() if v}


def reduce_native(row):
    return reduce_words((tuple(tuple(TO_LABEL[z] for z in t) for t in w), c) for w, c in row.items())


def marked(power, dotted):
    terms = []
    commutator = [((0, 1), 1), ((1, 0), -1)]
    for (a, u), (b, v) in product(commutator, repeat=2):
        terms.append((((3, 4),) * power + ((dotted + a + b),), u * v))
    return reduce_words(terms)


def calculate(colours, workers):
    weight = (3, 3, 2, 2, 0) if colours == 7 else (2, 2, 2, 2, 2)
    blocks = native.words_by_shape(10, weight)
    stable = [row for words in blocks.values() for row in native.singlets(words)]
    reduced = [reduce_native(r) for r in stable]
    tests = [marked(3, ()), marked(2, (3, 4))] if colours == 7 else []
    full_dimension = (dimension("SO7", weight, 10, 6) if colours == 7
                      else len(stable) - physical_image_relations(4, 2, colours))
    elementary = np.eye(colours**2, dtype=np.int64).reshape((colours,) * 4)
    basis = np.array([elementary[i, j] - elementary[j, i] for i in range(colours) for j in range(i)])
    evaluated, frame = integral_rows(stable, "SO" + str(colours), 2 * weight[0],
                                    full_dimension, 67291, workers, basis=basis)
    # A maximal-rank integer evaluation is injective on the finite-colour singlets.
    assert evaluated.rank() == full_dimension
    relation_columns, nullity = evaluated.transpose().nullspace()
    assert len(stable) - full_dimension == nullity
    relations = []
    for j in range(nullity):
        column = [int(relation_columns[i, j]) for i in range(len(stable))]
        divisor = gcd(*column)
        relations.append([v // divisor for v in column])
    kernel = fmpz_mat(relations) if relations else fmpz_mat(0, len(stable))
    assert kernel.rank() == nullity and (kernel * evaluated).is_zero()
    keys = sorted(set().union(*(v.keys() for v in reduced + tests)))
    formal = fmpz_mat([[r.get(k, 0) for k in keys] for r in reduced])
    projected = kernel * formal
    commuting = sorted({word for key in keys for sign, word in [abelian(key)] if sign})
    commute = fmpz_mat([[sign if word == col else 0 for col in commuting]
                        for key in keys for sign, word in [abelian(key)]])
    partial = formal.rank() - projected.rank()
    full = (formal * commute).rank() - (projected * commute).rank()
    assert 0 <= full <= partial
    directions = [[r.get(k, 0) for k in keys] for r in tests]
    original = projected.rank()
    def augmented(test):
        return fmpz_mat(projected.tolist() + test).rank()
    ranks = [augmented([r]) for r in directions]
    both = augmented(directions) if directions else original
    # Any dual vector returned here annihilates every specialized trace identity.
    annihilator, n = projected.nullspace()
    values = fmpz_mat(directions) * annihilator if directions else None
    witnesses = []
    for j in range(n if directions else 0):
        left, right = int(values[0, j]), int(values[1, j])
        if left or right:
            v = [int(annihilator[i, j]) for i in range(len(keys))]
            divisor = gcd(*v)
            witnesses.append(dict(pairings=[left // divisor, right // divisor],
                                  coordinates=[int(x // divisor) for x in v]))
    rref, rank = fmpq_mat(projected).rref()
    return dict(group=f"SO({colours})", weight=weight, stable_singlets=len(stable),
                finite_singlets=full_dimension, finite_relations=nullity,
                reduced_monomials=len(keys), reduced_singlet_rank=formal.rank(),
                projected_relation_rank=original, augmented_ranks=ranks, joint_rank=both,
                native_marked_test_rank=both - original, frame=frame,
                ghost_quotient_dimension=partial, all_commuting_quotient_dimension=full,
                entire_scalar_commutator_image_dimension=partial - full,
                reduced_words=keys, independent_relation_rows=[[str(x) for x in r] for r in rref.tolist()[:rank]], marked_rows=directions,
                dual_witnesses=witnesses)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--colours", type=int, choices=(5, 7), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("assertions and a positive worker count are required")
    result = calculate(args.colours, args.workers)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print({k: result[k] for k in ("stable_singlets", "finite_singlets", "finite_relations",
                                  "projected_relation_rank", "augmented_ranks", "joint_rank",
                                  "entire_scalar_commutator_image_dimension")})
