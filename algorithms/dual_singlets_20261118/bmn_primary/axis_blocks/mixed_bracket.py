"""First differential of the auxiliary axis filtration on two mixed quadratic jets."""

from collections import Counter, defaultdict, deque
from functools import cache
from itertools import combinations, combinations_with_replacement

from flint import fmpz_mat

from g2_terminal import cross_image, one_axis_cycles
from transverse import colour_data
from bmn_primary_index import G2_WEYL, g2_dominant


HIGHEST = ((4, 2), (6, 0), (8, 0))


@cache
def columns_by_weight():
    cycles = one_axis_cycles()
    weights, _ = colour_data("G2")
    wanted = {g2_dominant((x+1-a, y+1-b)) for x, y in HIGHEST for a, b in G2_WEYL}
    spectators = defaultdict(list)
    for a, b in combinations(range(len(weights)), 2):
        spectators[weights[a][0]+weights[b][0], weights[a][1]+weights[b][1]].append((a, b))
    result = defaultdict(list)
    for i, j in combinations_with_replacement(range(len(cycles)), 2):
        x, y = (a+b for a, b in zip(cycles[i][0], cycles[j][0]))
        for wx, wy in wanted:
            for a, b in spectators.get((wx-x, wy-y), ()):
                result[wx, wy].append((i, j, a, b))
    return result


def differential(i, j, a, b):
    result = {}
    for (left, right, z), value in cross_image(i, j).items():
        if z not in (a, b):
            key = (left, right, tuple(sorted((a, b, z))))
            # Once a,b are fixed the third colour z is uniquely recovered.
            result[key] = (-1) ** ((a < z) + (b < z)) * value
    return result


def cubic_differential(i, j, a, b):
    """Merge one first jet on the first axis with either quadratic mixed jet."""
    cycles = one_axis_cycles()
    _, bracket = colour_data("G2")
    result = Counter()
    pairs = [(i, j)] if i == j else [(i, j), (j, i)]
    for left, right in pairs:
        for u, n in cycles[left][1]:
            for v, m in cycles[right][1]:
                for index, x in enumerate(u):
                    rest = u[:index] + u[index+1:]
                    for slot, (old, kept) in enumerate(((a, b), (b, a))):
                        for z, coefficient in bracket[x, old]:
                            result[rest, v, kept, z] += (-1) ** (index+slot) * n*m*coefficient
    return {w: n for w, n in result.items() if n}


def sparse_rank(columns):
    """Remove forced pivots before exact integer rank on the remaining rows."""
    rows = defaultdict(dict)
    for j, column in enumerate(columns):
        for w, coefficient in column.items():
            rows[w][j] = coefficient
    row_count = len(rows)
    pending = deque(w for w, r in rows.items() if len(r) == 1)
    pivots = set()
    while pending:
        w = pending.popleft()
        if len(rows[w]) != 1:
            continue
        j = next(iter(rows[w]))
        pivots.add(j)
        for key in columns[j]:
            del rows[key][j]
            if len(rows[key]) == 1:
                pending.append(key)
    remaining = [r for r in rows.values() if r]
    indices = {j: k for k, j in enumerate(sorted({j for r in remaining for j in r}))}
    matrix = fmpz_mat(len(remaining), len(indices))
    for i, row in enumerate(remaining):
        for j, value in row.items():
            matrix[i, indices[j]] = value
    return row_count, len(pivots) + matrix.rank()


@cache
def rank_at(weight, maximum_cost=4):
    cols = columns_by_weight().get(weight, ())
    images = [differential(*c) for c in cols]
    if maximum_cost >= 5:
        for col, arguments in zip(images, cols, strict=True):
            col.update(cubic_differential(*arguments))
    if maximum_cost >= 8:
        _, bracket = colour_data("G2")
        for col, (i, j, a, b) in zip(images, cols, strict=True):
            col.update({("quartic", i, j, z): n for z, n in bracket[a, b]})
    sites, rank = sparse_rank(images)
    return [*weight, len(cols), sites, rank, len(cols) - rank]


def calculate(maximum_cost=4):
    assert maximum_cost in (4, 5, 8)
    records = [rank_at(w, maximum_cost) for w in sorted(columns_by_weight())]
    multiplicities = []
    for x, y in HIGHEST:
        source = sum(s*rank_at(g2_dominant((x+1-a, y+1-b)), maximum_cost)[2] for (a, b), s in G2_WEYL.items())
        survivors = sum(s*rank_at(g2_dominant((x+1-a, y+1-b)), maximum_cost)[5] for (a, b), s in G2_WEYL.items())
        assert 0 <= survivors <= source
        multiplicities.append([x, y, source, survivors])
    return dict(weight_ranks=records, highest_weight_kernels=multiplicities)
