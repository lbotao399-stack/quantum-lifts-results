"""Singlets with a single contracted pair of Albert quartics and ten adjoint inputs."""
from collections import Counter
from functools import cache
from itertools import combinations, permutations, product
from math import gcd

from trace_ideal import echelon


# The scalar letters are 0,1,2 and the odd dotted letters are 3,4.
# The first three templates arise by splitting a nine-input covariant.
TEMPLATES = (
    ("nested", (((0, 1), 2), 3, 4, 5), (6, 7, 8, 9), ()),
    ("same_arm", ((0, 1), (2, 3), 4, 5), (6, 7, 8, 9), ()),
    ("opposite_arms", ((0, 1), 2, 3, 4), ((5, 6), 7, 8, 9), ()),
    ("quadratic_spectator", (0, 1, 2, 3), (4, 5, 6, 7), (8, 9)),
)

LETTERS = ((0, 0, (0,)), (0, 0, (1,)), (0, 0, (2,)), (1, 0, ()), (0, 1, ()),
           (2, 0, ()), (1, 1, ()), (0, 2, ()),
           *((a, b, (i,)) for a, b in ((1, 0), (0, 1)) for i in range(3)))


def sign(sequence):
    return (-1) ** sum(x > y for i, x in enumerate(sequence) for y in sequence[i + 1:])


@cache
def degree(node):
    if node[0] == 0:
        return int(len(LETTERS[node[1]][2]) % 2 == 0)
    return sum(degree(x) for x in node[1:])


def ordered(nodes, graded, antisymmetric=False):
    entries = tuple(nodes)
    odd = tuple(node for node in entries if degree(node) % 2)
    if not graded:
        return 1, tuple(sorted(entries))
    if antisymmetric:
        left, right = entries
        if left == right and not degree(left) % 2:
            return 0, entries
        phase = -1 if left > right else 1
    else:
        if len(set(odd)) != len(odd):
            return 0, entries
        phase = 1
    return phase * sign(odd), tuple(sorted(entries))


@cache
def primitive(a, b):
    from fortuity.full.letters import Q_letter
    x, y = LETTERS[a], LETTERS[b]
    z = (x[0]+y[0], x[1]+y[1], tuple(sorted(x[2]+y[2])))
    terms = {(u, v): c for c, u, v in Q_letter(z)}
    coefficient = (2 if a == b else 1) * terms[x, y]
    assert {k: coefficient*v for k, v in terms.items()} == (
        {(x, x): 2} if x == y else {(x,y): 1, (y,x): -(-1)**(degree((0,a))*degree((0,b)))})
    return coefficient, (0, LETTERS.index(z))


def colour_graph(template, labels, graded=True, marked=None):
    def build(tree):
        if isinstance(tree, int):
            return 1, (0, labels[tree])
        if tree == marked:
            if not graded:
                pair = sorted((0, labels[i]) for i in tree)
                return 1, (3, *pair)
            sign_before = (-1) ** sum(x >= 3 for x in labels[:tree[0]])
            coefficient, node = primitive(labels[tree[0]], labels[tree[1]])
            return sign_before * coefficient, node
        a, u = build(tree[0])
        b, v = build(tree[1])
        c, children = ordered((u, v), graded, antisymmetric=True)
        return a * b * c, (1, *children)

    def arm(sites):
        nodes = [build(x) for x in sites]
        phase, children = ordered((v for _, v in nodes), graded)
        for s, _ in nodes:
            phase *= s
        return phase, (2, *children)

    la, a = arm(template[1])
    lb, b = arm(template[2])
    lc, arms = ordered((a, b), graded)
    ls, spectator = arm(template[3])
    return la * lb * lc * ls, (*arms, spectator)


def innermost_pairs(tree):
    if isinstance(tree, int):
        return
    if all(isinstance(x, int) for x in tree):
        yield tree
    else:
        for child in tree:
            yield from innermost_pairs(child)


@cache
def assignments(dotted, spin_pair, flavour_triple):
    scalar = sorted(set(range(10)) - set(dotted))
    spin = (spin_pair, sorted(set(dotted) - set(spin_pair)))
    flavour = (flavour_triple, sorted(set(scalar) - set(flavour_triple)))
    pairs = list(permutations((3, 4)))
    triples = list(permutations(range(3)))
    out = []
    for a, b, c, d in product(pairs, pairs, triples, triples):
        values = [0] * 10
        for index, label in zip((*spin[0], *spin[1], *flavour[0], *flavour[1]), (*a, *b, *c, *d)):
            values[index] = label
        out.append((tuple(values), sign(a) * sign(b) * sign(c) * sign(d)))
    return out


def singlets(sources=False):
    candidates = []
    for template in TEMPLATES:
        innermost = [pair for arm in template[1:3] for tree in arm for pair in innermost_pairs(tree)]
        layouts = set()
        for marked, dotted in product(innermost if sources else (None,), combinations(range(10), 4)):
            if marked and not set(marked).intersection(dotted):
                continue
            symbol = tuple(3 if j in dotted else 0 for j in range(10))
            _, layout = colour_graph(template, symbol, graded=False, marked=marked)
            if layout in layouts:
                continue
            layouts.add(layout)
            inactive = any(x in dotted or y in dotted for x, y in innermost)
            scalar = sorted(set(range(10)) - set(dotted))
            for partner in dotted[1:]:
                for tail in combinations(scalar[1:], 2):
                    row = Counter()
                    pair, triple = (dotted[0], partner), (scalar[0], *tail)
                    for labels, coefficient in assignments(dotted, pair, triple):
                        phase, graph = colour_graph(template, labels, marked=marked)
                        if phase:
                            row[graph] += phase * coefficient
                    row = {g: v for g, v in row.items() if v}
                    if row:
                        divisor = gcd(*row.values())
                        data = (template[0], dotted, pair, triple, divisor)
                        if marked:
                            data += (marked,)
                        candidates.append((inactive, data, {g:v // divisor for g, v in row.items()}))
    basis, rows = {}, []
    for inactive, data, row in sorted(candidates, key=lambda item: not item[0]):
        previous = len(basis)
        echelon(basis, row)
        if len(basis) > previous:
            rows.append((inactive, data, row))
    return rows


if __name__ == "__main__":
    from collections import Counter
    rows = singlets()
    print(Counter((inactive, recipe[0]) for inactive, recipe, _ in rows))
    print(len({g for _, _, row in rows for g in row}))
