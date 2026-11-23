"""Compute the first possible transverse source symbol at commuting F4 root quadruples."""

from collections import Counter, defaultdict
from itertools import combinations
import argparse
import json
from pathlib import Path

from flint import fmpz_mat

from colour import Colour
from commuting_source import simple_coordinates
from nilpotent_flag_orders import quad_orbits


def add(*vectors):
    return tuple(map(sum, zip(*vectors)))


def neg(vector):
    return tuple(-x for x in vector)


def weight_blocks(colour, roots, quartet):
    spaces = [defaultdict(list) for _ in range(5)]
    for degree, space in enumerate(spaces):
        for axes in combinations(range(4), degree):
            shift = tuple(sum(roots[quartet[i]][j] for i in axes) for j in range(4))
            for index, weight in enumerate(roots):
                space[add(weight, neg(shift))].append((axes, index))
    for weight in sorted(set().union(*spaces)):
        dims = [len(space.get(weight, ())) for space in spaces]
        matrices = []
        for degree in range(4):
            target = spaces[degree + 1].get(weight, ())
            rows = {key: i for i, key in enumerate(target)}
            block = fmpz_mat(dims[degree + 1], dims[degree])
            for j, (axes, source) in enumerate(spaces[degree].get(weight, ())):
                for i, anchor in enumerate(quartet):
                    if i in axes:
                        continue
                    ordered = tuple(sorted((*axes, i)))
                    sign = (-1) ** ordered.index(i)
                    for k, coeff in colour.br[anchor][source].items():
                        block[rows[ordered, k], j] = sign * coeff
            if matrices:
                assert (block * matrices[-1]).is_zero()
            matrices.append(block)
        yield weight, [space.get(weight, ()) for space in spaces], matrices


def tangent(colour, roots, quartet):
    cohomology = [Counter() for _ in range(5)]
    ranks = [0] * 4
    for weight, bases, matrices in weight_blocks(colour, roots, quartet):
        dims = list(map(len, bases))
        local = [b.rank() for b in matrices]
        for i, rank in enumerate(local):
            ranks[i] += rank
        for i, (before, after) in enumerate(zip([0] + local, local + [0])):
            multiplicity = dims[i] - before - after
            assert multiplicity >= 0
            if multiplicity:
                cohomology[i][weight] = multiplicity
    volume = add(*(roots[i] for i in quartet))
    for i in range(5):
        dual = {neg(add(volume, w)): n for w, n in cohomology[i].items()}
        assert dual == cohomology[4 - i]
    return cohomology, ranks


def encode(w):
    # All sums queried below have at most 25 adjoint-root contributions per axis.
    # Use exact addition and an injective digit range during each search.
    return sum(x * 2048**i for i, x in enumerate(w))


def first_orders(first, second, targets, limit):
    tangent_weights = {encode(w): w for w in first}
    mark_weights = {encode(w): w for w in second}
    max_digit = max(abs(x) for w in (*first, *second, *targets) for x in w)
    assert (limit + 3) * max_digit < 1024
    half = [{0}]
    ordered_half = [(0,)]

    def layer(k):
        while len(half) <= k:
            half.append({a + b for a in half[-1] for b in tangent_weights})
        return half[k]

    results = {}
    for target in targets:
        wanted = encode(neg(target))
        for k in range(limit + 1):
            right = layer(k - k // 2)
            while len(ordered_half) <= k // 2:
                ordered_half.append(tuple(sorted(layer(len(ordered_half)))))
            witness = next(
                (
                    (m, x, wanted - m - x)
                    for m in mark_weights
                    for x in ordered_half[k // 2]
                    if wanted - m - x in right
                ),
                None,
            )
            if witness is not None:
                m, x, y = witness

                def unpack(value, size):
                    out = []
                    for j in range(size, 0, -1):
                        v = next(
                            v for v in tangent_weights if value - v in layer(j - 1)
                        )
                        out.append(tangent_weights[v])
                        value -= v
                    assert value == 0
                    return out

                directions = unpack(x, k // 2) + unpack(y, k - k // 2)
                assert not any(add(target, mark_weights[m], *directions))
                results[target] = dict(
                    order=k, mark=mark_weights[m], tangent=directions
                )
                break
        else:
            results[target] = dict(order=None, excluded_through=limit)
    return results


def calculate(limit):
    c = Colour()
    roots = simple_coordinates(c)
    raw = json.loads(
        Path(__file__).with_name("nilpotent_flag_orders.json").read_bytes()
    )
    previous = {
        tuple(row["root_colours"]): row["first_allowed_order"] for row in raw["orbits"]
    }
    result = []
    for quartet, _, _ in quad_orbits(c, roots):
        assert all(not c.br[i][j] for i, j in combinations(quartet, 2))
        cohomology, ranks = tangent(c, roots, quartet)
        total = add(*(roots[i] for i in quartet))
        targets = [
            tuple(4 * s - 2 * t for s, t in zip(total, roots[i])) for i in quartet
        ]
        orders = first_orders(cohomology[1], cohomology[2], targets, limit)
        tests = [orders[target] for target in targets]
        assert all(t["order"] is None or t["order"] >= previous[quartet] for t in tests)
        result.append(
            dict(
                root_colours=quartet,
                ranks=ranks,
                tangent_cohomology_dimensions=[sum(h.values()) for h in cohomology],
                linear_bound_before=previous[quartet],
                short_axis_tests=tests,
            )
        )
    return dict(
        transverse_cutoff=limit,
        group="F4",
        ordinary_shape=(4, 4, 4, 2),
        roots="native simple-root coordinates",
        orbits=result,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or not 0 <= args.order <= 12:
        parser.error("enable assertions and specify an order between zero and twelve")
    data = calculate(args.order)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
    print(
        Counter(
            min(
                (t["order"] for t in r["short_axis_tests"] if t["order"] is not None),
                default=None,
            )
            for r in data["orbits"]
        )
    )
