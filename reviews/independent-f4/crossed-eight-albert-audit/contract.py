"""Contract an Albert graph over Q, with a fixed exterior order on each edge list."""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import lcm
import numpy as np


def parity(left, right):
    inversions = sum((right & ((1 << k) - 1)).bit_count()
                     for k in range(left.bit_length()) if left & (1 << k))
    return (-1) ** inversions


def absolute_bound(tensor):
    return max(abs(int(tensor.min())), abs(int(tensor.max())))


def matrix_product(a, b):
    assert a.shape[-1] * absolute_bound(a) * absolute_bound(b) < 2**63
    return a @ b


def form(n, entries):
    divisor = lcm(*(Fraction(v).denominator for _, _, v in entries))
    a = np.zeros((n, n), dtype=np.int64)
    for i, j, v in entries:
        a[i, j] = int(Fraction(v) * divisor)
    return divisor, a


class Colour:
    def __init__(self, data):
        one, h = form(26, data['representationInvariantForm'])
        unit, b = form(52, data['metric'])
        self.d, inverse = form(52, data['inverseMetric'])
        assert one == unit == 1
        hi = 3*h.copy()
        hi[:2, :2] = [[2, 1], [1, 2]]
        assert np.array_equal(h @ hi, 3*np.eye(26, dtype=np.int64))
        assert np.array_equal(b @ inverse, self.d*np.eye(52, dtype=np.int64))
        self.propagator = {26: hi, 52: inverse}
        self.cubic = np.zeros((26, 26, 26), dtype=np.int64)
        for i, j, k, value in data['cubic']:
            self.cubic[i, j, k] = value
        rep = np.zeros((52, 26, 26), dtype=np.int64)
        for a, sparse in enumerate(data['representation26']):
            for i, j, value in sparse:
                rep[a, i, j] = int(value)
        adj = np.zeros((52, 52, 52), dtype=np.int64)
        for a, j, row in data['bracketPairs']:
            for i, value in row:
                adj[a, i, j], adj[j, i, a] = int(value), -int(value)
        self.vertex = {26: matrix_product(rep.transpose(0, 2, 1), h),
                       52: matrix_product(adj.transpose(0, 2, 1), b)}
        self.generator = {}
        for n, matrices in ((26, rep), (52, adj)):
            raised = matrix_product(inverse, matrices.reshape(52, -1))
            self.generator[n] = raised.reshape(52, n, n).transpose(0, 2, 1)


def act(tensor, matrix, axis):
    aligned = np.moveaxis(tensor, axis, -1)
    result = matrix_product(aligned.reshape(-1, aligned.shape[-1]), matrix)
    return np.moveaxis(result.reshape(aligned.shape), -1, axis)


def contract(colour, graph, bits, tree, poles, canonical=True):
    edges, dims = graph['E'], graph['d']
    plus, minus = graph['cp'], graph['cm']
    occupied = sorted(set(plus) | set(minus))
    owner = {e: edges[e][s] for e, s in zip(occupied, bits, strict=True)}
    kernels = {}
    for e in occupied:
        n = dims[e]
        for i, j in product(range(6) if e in plus else [-1],
                            range(6) if e in minus else [-1]):
            h = colour.propagator[n]
            if i >= 0:
                h = matrix_product(h, colour.generator[n][poles[0][i]])
            if j >= 0:
                h = matrix_product(h, colour.generator[n][poles[1][j]])
            kernels[e, i, j] = h if owner[e] == edges[e][0] else h.T

    def leaf(vertex):
        incident = tuple(e for e, ends in enumerate(edges) if vertex in ends)
        if vertex < graph['cubics']:
            tensor = colour.cubic
        else:
            rail = [e for e in incident if dims[e] == 26]
            fan = next(e for e in incident if e not in rail) if rail else len(edges) - 1
            arms = [e for e in incident if e != fan]
            entering = next(e for e in arms if edges[e][1] == vertex)
            exiting = next(e for e in arms if edges[e][0] == vertex)
            slots = [fan, entering, exiting]
            tensor = colour.vertex[dims[entering]].transpose([slots.index(e) for e in incident])
        terms = {(0, 0): tensor}
        source = [0, 0]
        for axis, e in enumerate(incident):
            if e not in owner:
                if edges[e][0] == vertex:
                    terms = {m: act(t, colour.propagator[dims[e]], axis) for m, t in terms.items()}
                continue
            if owner[e] != vertex:
                continue
            slots = [1 << group.index(e) if e in group else 0 for group in (plus, minus)]
            edge_sign = parity(source[0], slots[0]) * parity(source[1], slots[1]) if canonical else 1
            source = [a | b for a, b in zip(source, slots)]
            updated = {}
            for (u, v), t in terms.items():
                for i, j in product(range(6) if slots[0] else [-1], range(6) if slots[1] else [-1]):
                    x, y = (1 << i if i >= 0 else 0), (1 << j if j >= 0 else 0)
                    if u & x or v & y:
                        continue
                    k = u | x, v | y
                    val = edge_sign * parity(u, x) * parity(v, y) * act(t, kernels[e, i, j], axis)
                    if k in updated:
                        assert absolute_bound(updated[k]) + absolute_bound(val) < 2**63
                        updated[k] += val
                    else:
                        updated[k] = val
            terms = updated
        sparse = {}
        for mask, t in terms.items():
            index = np.nonzero(t)
            if len(index[0]):
                sparse[mask] = dict(zip(zip(*index), map(int, t[index])))
        return incident, source, sparse

    def group(values, axes, remaining, shared):
        outside = [axes.index(e) for e in remaining]
        inside = [axes.index(e) for e in shared]
        grouped = {}
        for mask, row in values.items():
            block = defaultdict(list)
            for ix, val in row.items():
                key = tuple(ix[k] for k in inside)
                block[key].append((tuple(ix[k] for k in outside), val))
            grouped[mask] = block
        return grouped

    def visit(node):
        if isinstance(node, int):
            return leaf(node)
        a, sa, x = visit(node[0])
        b, sb, y = visit(node[1])
        shared = [e for e in a if e in b]
        l, r = [e for e in a if e not in shared], [e for e in b if e not in shared]
        left, right = group(x, a, l, shared), group(y, b, r, shared)
        edge_sign = parity(sa[0], sb[0]) * parity(sa[1], sb[1]) if canonical else 1
        out = defaultdict(lambda: defaultdict(int))
        for (u, v), first in left.items():
            for (i, j), second in right.items():
                if u & i or v & j:
                    continue
                s = edge_sign * parity(u, i) * parity(v, j)
                dest = out[u | i, v | j]
                for k in (first if len(first) <= len(second) else second):
                    if k in first and k in second:
                        for ix, value in first[k]:
                            value *= s
                            for iy, other in second[k]:
                                dest[ix + iy] += value * other
        clean = {k: {i: val for i, val in vals.items() if val} for k, vals in out.items()}
        return tuple(l + r), [sa[0] | sb[0], sa[1] | sb[1]], clean

    axes, source, values = visit(tree)
    assert not axes and source == [63, 63]
    denominator = 3**dims.count(26) * colour.d**(dims.count(52) + 12)
    return Fraction(values.get((63, 63), {}).get((), 0), denominator)
