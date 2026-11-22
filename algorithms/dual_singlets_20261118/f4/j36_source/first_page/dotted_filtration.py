"""Classify one-scalar dotted cycles by testing weights forced by the scalar-jet filtration."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, groupby
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from f4_singlet_characters import CARTAN, ZERO, multiply
from dotted_syzygies import Dotted
from mixed_syzygies import PRIME

SIMPLE = CARTAN["F4"]
# Twice the Gram matrix of fundamental weights, with long roots of square two.
GRAM2 = ((4, 6, 4, 2), (6, 12, 8, 4), (4, 8, 6, 3), (2, 4, 3, 2))


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def difference(x, y):
    return tuple(a-b for a, b in zip(x, y))


def root_coordinates(labels):
    lengths2 = (2, 2, 1, 1)
    coords = tuple(Fraction(dot(labels, col), length)
                   for col, length in zip(zip(*GRAM2), lengths2))
    assert all(v.denominator == 1 for v in coords)
    return tuple(int(v) for v in coords)


def casimir2(labels):
    return sum(x*dot(col, tuple(y+2 for y in labels))
               for x, col in zip(labels, GRAM2))


@cache
def exterior(k):
    """Newton identities in the integral Weyl-character ring."""
    if k == 0:
        return {ZERO: 1}
    result = Counter()
    for j in range(1, k+1):
        for highest, coeff in exterior(k-j).items():
            for target, mult in multiply("F4", highest, j):
                result[target] += (-1)**(j-1)*coeff*mult
    assert all(n >= 0 and n % k == 0 for n in result.values())
    return {w: v//k for w, v in result.items() if v}


def tensor_adjoint(character):
    result = Counter()
    for highest, n in character.items():
        for target, mult in multiply("F4", highest, 1):
            result[target] += n*mult
    return {w: v for w, v in result.items() if v}


def harmonic(character, t_degree):
    # The positive-current Laplacian on an irreducible is proportional to
    # 2 h^vee times the loop degree minus its finite quadratic Casimir.
    return {w: v for w, v in character.items() if casimir2(w) == 36*t_degree}


def pure_current(t_degree, excess):
    if not excess:
        return harmonic(exterior(t_degree), t_degree)
    assert excess == 1 and t_degree >= 2
    return harmonic(tensor_adjoint(exterior(t_degree-2)), t_degree)


@cache
def lowering_path(highest, target):
    if highest == target:
        return ()
    delta = root_coordinates(difference(highest, target))
    if min(delta) >= 0:
        for i, room in enumerate(delta):
            if room and highest[i] > 0:
                path = lowering_path(difference(highest, SIMPLE[i]), target)
                if path is not None:
                    return (i, *path)
    return None


class Filtered(Dotted):
    @cache
    def colour_groups(self, size):
        by_weight = defaultdict(list)
        for indices in combinations(range(52), size):
            weight = tuple(map(sum, zip(*(self.c.weights[i] for i in indices))))
            by_weight[weight].append(indices)
        return by_weight

    def assignments(self, slots, weight):
        # Keep equal-mode odd colours strictly ordered. No table exceeds C(52,3).
        chunks = []
        for mode, it in groupby(slots):
            count = len(list(it))
            while count:
                k = min(3, count)
                chunks.append((mode, k))
                count -= k

        def choose(pos, remaining, prefix):
            if pos == len(chunks):
                if remaining == ZERO:
                    yield prefix
                return
            mode, size = chunks[pos]
            groups = self.colour_groups(size)
            same_mode = pos > 0 and chunks[pos-1][0] == mode
            choices = (remaining,) if pos+1 == len(chunks) else groups
            for w in choices:
                rest = difference(remaining, w)
                if pos+2 == len(chunks) and rest not in self.colour_groups(chunks[-1][1]):
                    continue
                for colours in groups.get(w, ()):
                    if not same_mode or prefix[-1] < colours[0]:
                        yield from choose(pos+1, rest, prefix+colours)

        yield from choose(0, weight, ())

    def native_weight(self, labels):
        simple = [self.c.weights[i] for i in self.c.data["raising"]]
        coords = root_coordinates(labels)
        w = tuple(dot(coords, col) for col in zip(*simple))
        for label, e, f in zip(labels, self.c.data["raising"], self.c.data["lowering"]):
            h = self.c.br[e][f]
            assert 2*sum(w[i]*v for i, v in h.items()) == label*sum(
                self.c.weights[e][i]*v for i, v in h.items())
        return w


def packed(character):
    return [[list(w), n] for w, n in sorted(character.items())]


def calculate(n, perform_ranks=True):
    engine = Filtered()
    for i, row in enumerate(SIMPLE):
        assert root_coordinates(row) == tuple(int(i == j) for j in range(4))
    # A scalar with zero derivatives leaves one excess in the pure ghost
    # complex. A scalar with one derivative leaves only first ghosts.
    ghost_excess = pure_current(n-1, 1)
    first_ghosts = pure_current(n-2, 0)
    candidates = sorted(tensor_adjoint(Counter(ghost_excess)+Counter(first_ghosts)))
    minima = [w for w in candidates if all(
        w == v or lowering_path(w, v) is None for v in candidates)]
    paths = []
    for w in candidates:
        bound = next(v for v in minima if lowering_path(w, v) is not None)
        steps = lowering_path(w, bound)
        cur = w
        for i in steps:
            assert cur[i] > 0
            cur = difference(cur, SIMPLE[i])
        assert cur == bound
        paths.append(dict(highest=w, tested=bound, simple_lowerings=steps))
    weights = []
    for w in minima:
        native = engine.native_weight(w)
        row = dict(labels=w, native=native)
        if perform_ranks:
            row.update(engine.slice(n, native))
        weights.append(row)
    return dict(slots=n, scalar_degree=1, dotted_degree=n-1,
                all_colour_columns=52*comb(52, n-2)+52**2*comb(52, n-3),
                pure_ghost_one_excess=packed(ghost_excess), first_ghosts=packed(first_ghosts),
                highest_weight_paths=paths, prime=PRIME, tested_weights=weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("slots", type=int, choices=range(4, 10))
    parser.add_argument("--support-only", action="store_true")
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    print(json.dumps(calculate(args.slots, not args.support_only), indent=2))
