"""Rational colour and permutation characters for a block of at most four Taylor slots."""

from collections import Counter, defaultdict
from functools import cache
from itertools import combinations, permutations, product
from pathlib import Path
import sys

from flint import fmpz_mat

PRIMARY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PRIMARY))
sys.path.insert(0, str(PRIMARY.parent))
from b2_schur_complex import character, partitions
from bmn_primary_index import g2_orbit
from symplectic_axis_bounds import structure
from g2_colours import structure as g2_structure


@cache
def colour_data(group="USp(4)"):
    return {"USp(4)": structure, "G2": g2_structure}[group]()


def dominant(x, y, group):
    return x >= y >= 0 if group == "USp(4)" else min(x, y) >= 0


@cache
def weight_tables(weight, group="USp(4)"):
    """CE homology at an ordinary Taylor weight, resolved by dominant torus weights."""
    colours, brackets = colour_data(group)
    modes = [v for v in product(*(range(n + 1) for n in weight)) if any(v)]
    mode_index = {v: i for i, v in enumerate(modes)}
    spaces = defaultdict(list)
    ncol = len(colours)

    def build(left, first, word, x, y):
        if not any(left):
            if dominant(x, y, group):
                spaces[x, y, len(word)].append(word)
            return
        for mode in range(first // ncol, len(modes)):
            residual = tuple(a - b for a, b in zip(left, modes[mode]))
            if min(residual) < 0:
                continue
            for colour in range(max(0, first - mode * ncol), ncol):
                letter = ncol * mode + colour
                dx, dy = colours[colour]
                build(residual, letter + 1, word + (letter,), x + dx, y + dy)

    build(weight, 0, (), 0, 0)

    @cache
    def join(a, b):
        mode = tuple(i + j for i, j in zip(modes[a // ncol], modes[b // ncol]))
        return [(ncol * mode_index[mode] + c, n) for c, n in brackets[a % ncol, b % ncol]]

    cochains, result = Counter(), Counter()
    for x, y in sorted({(a, b) for a, b, _ in spaces}):
        previous = fmpz_mat(0, 0)
        previous_rank = 0
        for length in range(sum(weight) + 2):
            source = spaces.get((x, y, length), ())
            target = spaces.get((x, y, length - 1), ())
            indices = {v: k for k, v in enumerate(target)}
            boundary = fmpz_mat(len(target), len(source))
            for col, word in enumerate(source):
                for i, j in combinations(range(length), 2):
                    rest = word[:i] + word[i + 1:j] + word[j + 1:]
                    for letter, value in join(word[i], word[j]):
                        if letter not in rest:
                            row = indices[tuple(sorted((*rest, letter)))]
                            sign = (-1) ** (i + j - 1 + sum(a < letter for a in rest))
                            boundary[row, col] += sign * value
            assert (previous * boundary).is_zero()
            rank = boundary.rank()
            dimension = len(target) - rank - previous_rank
            assert dimension >= 0
            if target:
                cochains[x, y, length - 1] = len(target)
            if dimension:
                result[x, y, length - 1] = dimension
            previous, previous_rank = boundary, rank
    return cochains, result


@cache
def kostka(shape, weight):
    """Remove the horizontal strip occupied by the largest tableau entry."""
    if not weight:
        return int(not shape)
    if not shape or sum(shape) != sum(weight):
        return 0
    strips = [range(b, a + 1) for a, b in zip(shape, (*shape[1:], 0))]
    return sum(kostka(tuple(v for v in smaller if v), weight[:-1])
               for smaller in product(*strips) if sum(smaller) == sum(shape) - weight[-1])


@cache
def schur_table(degree, homology, group="USp(4)", rows=None):
    shapes = sorted((p[::-1] for p in partitions(degree) if rows is None or len(p) <= rows), reverse=True)
    result = {}
    for shape in shapes:
        multiplicities = weight_tables(shape, group)[homology].copy()
        for earlier, terms in result.items():
            factor = kostka(earlier, shape)
            if factor:
                multiplicities.subtract({w: factor * v for w, v in terms.items()})
        assert all(v >= 0 for v in multiplicities.values()) and kostka(shape, shape) == 1
        result[shape] = +multiplicities
    return result


@cache
def signed_orbit(x, y, group="USp(4)"):
    if group == "G2":
        return tuple(g2_orbit((x, y)))
    return {(i * a, j * b) for a, b in permutations((x, y)) for i, j in product((-1, 1), repeat=2)}


@cache
def low_character(cycles, length, homology, group="USp(4)", rows=None):
    """Trace of the given permutation by colour weight on the multilinear homology."""
    if not cycles:
        return {(0, 0): 1} if length == 0 else {}
    terms = Counter()
    for shape, weights in schur_table(sum(cycles), homology, group, rows).items():
        factor = character(shape, cycles)
        for (x, y, degree), value in weights.items():
            if degree == length and factor:
                for w in signed_orbit(x, y, group):
                    terms[w] += factor * value
    return {w: n for w, n in terms.items() if n}
