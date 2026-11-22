"""Multilinear classical SU(2) cochains with one second jet, in Taylor degrees 3, 4, 5, 7."""
from functools import cache
from itertools import combinations, product
import json

import numpy as np
from flint import fmpz_mat


@cache
def matchings(slots):
    if not slots:
        return ((),)
    first, *rest = slots
    return tuple(((first, j),) + tail for j in rest
                 for tail in matchings(tuple(k for k in rest if k != j)))


def invariant_tensors(size):
    slots = tuple(range(size))
    if size % 2 == 0:
        return [((), pairs) for pairs in matchings(slots)]
    return [(triple, pairs) for triple in combinations(slots, 3)
            for pairs in matchings(tuple(i for i in slots if i not in triple))]


def orientation(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    return 1 if (j - i) % 3 == 1 else -1


def tensor_value(tensor, colours):
    triple, pairs = tensor
    if any(colours[i] != colours[j] for i, j in pairs):
        return 0
    return orientation(*(colours[i] for i in triple)) if triple else 1


def bracket(i, j):
    if i == j:
        return 0, 0
    k = 3 - i - j
    return k, orientation(i, j, k)


def with_bracket(tensor, colours, place, i, j):
    value, sign = bracket(i, j)
    if not sign:
        return 0
    entries = list(colours)
    entries[place] = value
    return sign * tensor_value(tensor, entries)


def probe(tensor, size):
    # Up to six vector slots: a separate axis for each metric picks just this tensor.
    triple, pairs = tensor
    cols = [0] * size
    if triple:
        assert size == 3
        cols[:] = range(3)
    for axis, (left, right) in enumerate(pairs):
        cols[left] = cols[right] = axis
    return cols


def integer_rank(matrix):
    if not matrix.size:
        return 0
    return fmpz_mat(matrix.tolist()).rank()


def run(n):
    labels = tuple(range(n))
    pairs = tuple(combinations(labels, 2))
    mid_tensors = invariant_tensors(n - 1)
    parent_tensors = invariant_tensors(n - 2)
    mid = [(p, t) for p in pairs for t in mid_tensors]
    singles = {p: [i for i in labels if i not in p] for p in pairs}
    probes = [probe(t, n - 1) for t in mid_tensors]
    assert [[tensor_value(t, v) for t in mid_tensors] for v in probes] == np.eye(
        len(mid_tensors), dtype=int).tolist()

    outputs = set()
    for colours in product(range(3), repeat=n):
        values = []
        for p in pairs:
            i, j = p
            k, sign = bracket(colours[i], colours[j])
            args = [colours[s] for s in singles[p]] + [k]
            values.extend(sign * tensor_value(t, args) for t in mid_tensors)
        if any(values):
            outputs.add(tuple(values))
    outgoing = np.array(sorted(outputs), dtype=np.int64).reshape(-1, len(mid))

    parents = []
    for triple in combinations(labels, 3):
        parents.extend(((triple,), t) for t in parent_tensors)
    for p, q in combinations(pairs, 2):
        if set(p).isdisjoint(q):
            parents.extend(((p, q), t) for t in parent_tensors)

    incoming = np.zeros((len(mid), len(parents)), dtype=np.int64)
    for i, (pair, t) in enumerate(mid):
        scalar = dict(zip(singles[pair], probe(t, n - 1)[:-1], strict=True))
        marked = probe(t, n - 1)[-1]
        for j, (blocks, tensor) in enumerate(parents):
            args = [scalar[s] for s in labels if all(s not in b for b in blocks)
                    and s in scalar]
            if len(blocks) == 1:
                triple = blocks[0]
                if not set(pair) <= set(triple):
                    continue
                at = next(k for k in triple if k not in pair)
                sign = (-1) ** triple.index(at)
                incoming[i, j] = sign * with_bracket(
                    tensor, args + [0], -1, scalar[at], marked)
            elif pair in blocks:
                other = blocks[1 - blocks.index(pair)]
                slot = blocks.index(other)
                incoming[i, j] = (-1) ** slot * with_bracket(
                    tensor, args + [marked, marked], len(args) + slot,
                    scalar[other[0]], scalar[other[1]])

    # Pair products times the two closed metric combinations for each labelled triple.
    column_ids = {slot: j for j, slot in enumerate(mid)}
    cycles = []
    if n % 2:
        for a, b, c in combinations(labels, 3):
            rest = tuple(i for i in labels if i not in (a, b, c))
            for spectator in matchings(rest):
                entries = []
                for scalar, pair in ((a, (b, c)), (b, (a, c)), (c, (a, b))):
                    positions = {k: j for j, k in enumerate(singles[pair])}
                    terms = [(positions[l], positions[r]) for l, r in spectator]
                    terms.append((positions[scalar], n - 2))
                    entries.append(column_ids[pair, ((), tuple(sorted(terms)))])
                for left, right in zip(entries, entries[1:]):
                    row = [0] * len(mid)
                    row[left] = row[right] = 1
                    cycles.append(row)
    cycles = np.array(cycles, dtype=np.int64).reshape(-1, len(mid)).T
    cartan = np.zeros((len(pairs), len(mid)), dtype=np.int64)
    for i, (p, t) in enumerate(mid):
        cartan[pairs.index(p), i] = tensor_value(t, [0] * (n - 1))
    assert not (outgoing @ incoming).any()
    assert not (outgoing @ cycles).any()
    assert not (cartan @ incoming).any()
    image = integer_rank(incoming)
    out = integer_rank(outgoing)
    visible = integer_rank(cartan @ cycles)
    cohomology = len(mid) - image - out
    return dict(degree=n, middle=len(mid), parent_columns=len(parents),
                incoming_rank=image, outgoing_rank=out, cohomology=cohomology,
                cartan_rank_of_cubic_cycles_and_products=visible,
                cartan_kernel_upper_bound=cohomology - visible)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    print(json.dumps([run(n) for n in (3, 4, 5, 7)], indent=2))
