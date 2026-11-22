"""The length-fifteen branch of the G2 BMN singlet at level 36."""

from collections import Counter
from itertools import product
import json

from paired_axes import (G2_WEYL, PRIMARY, axis, exterior, partner, projection, singlets,
                         tensor, trivial_multiplicity)
from g2_terminal import terminal_character
from mixed_bracket import calculate as mixed_ranks
from transverse import kostka


def occupations(weight, length):
    """Necessary pure-axis and mixed-jet lengths; a mixed jet uses at least two axes."""
    choices = [[(e, y) for e, y in axis("G2") if e <= x] for x in weight]
    rows = []
    for pure in product(*choices):
        residual = [v - e for v, (e, _) in zip(weight, pure)]
        mixed_length = length - sum(y for _, y in pure)
        if mixed_length < 0 or bool(any(residual)) != bool(mixed_length):
            continue
        # Each mixed jet needs at least one slot outside any single axis.
        if mixed_length <= min(sum(residual) // 2, sum(residual) - max(residual)):
            rows.append((pure, mixed_length))
    return rows


def first_two_pairs(mixed):
    traces = []
    high_energy, high_length = (3, 3) if mixed else (5, 4)
    for top, bottom in product((0, 1), repeat=2):
        first = ([("axis", 2, high_energy, high_length)] if top
                 else [("axis", 1, high_energy, high_length)] * 2)
        last = [("axis", 2, 3, 3)] if bottom else [("axis", 1, 3, 3)] * 2
        factors = first + last + ([("mixed", 1, 2, 0)] if mixed else [])
        value = (-1) ** (high_length * top + 3 * bottom) * singlets("G2", None, tuple(sorted(factors)))
        traces.append(value)
    a, b, c, d = traces
    assert (a - b - c + d) % 4 == 0
    return traces, (a - b - c + d) // 4


def highest(weights, x, y):
    return sum(s * weights.get((x + 1 - a, y + 1 - b), 0) for (a, b), s in G2_WEYL.items())


def second_bound():
    # Both mixed jets in the second possibility must use the two high axes.
    assert occupations((5, 5, 3, 3), 14) == [
        (((3, 3), (3, 3), (3, 3), (3, 3)), 2), (((5, 4), (5, 4), (3, 3), (3, 3)), 0)]
    assert occupations((4, 2), 6) == []
    traces = [first_two_pairs(mixed) for mixed in (False, True)]
    low, ranks = terminal_character()
    high_sectors = []
    for e, y, sign in ((5, 4, -1), (3, 3, 1)):
        ordinary = tensor("G2", None, (("axis", 1, e, y), ("axis", 1, e, y)))
        doubled = tensor("G2", None, (("axis", 2, e, y),))
        assert all((n + sign * doubled.get(w, 0)) % 2 == 0 for w, n in ordinary.items())
        square = {w: (n + sign * doubled.get(w, 0)) // 2 for w, n in ordinary.items()}
        multipliers = exterior("G2", 2) if e == 3 else {(0, 0): 1}
        weights = Counter()
        for (a, b), n in square.items():
            for (c, d), m in multipliers.items():
                weights[a+c, b+d] += n*m
        high_sectors.append(weights)
    channels = [[x, y, highest(low, x, y), *(highest(rep, x, y) for rep in high_sectors)]
                for x, y in sorted(low) if min(x, y) >= 0 and highest(low, x, y)]
    assert all(n > 0 and min(a, b) >= 0 for _, _, n, a, b in channels)
    # Weyl dimension formula in the short-root, long-root Dynkin label order.
    def dimension(x, y):
        return (x+1)*(y+1)*(x+y+2)*(x+2*y+3)*(x+3*y+4)*(2*x+3*y+5)//120
    assert sum(n * dimension(x, y) for x, y, n, _, _ in channels) == sum(low.values())
    bounds = [sum(n*(a+b) for _, _, n, a, b in channels)]
    mixed = {}
    for cost in (4, 5, 8):
        stage = mixed_ranks(cost)
        kernels = {(x, y): (before, after) for x, y, before, after in stage["highest_weight_kernels"]}
        bound = 0
        for x, y, n, plain, initial in channels:
            before, surviving = kernels.get((x, y), (0, 0))
            assert initial == before
            bound += n * (plain + surviving)
        assert bound <= bounds[-1]
        bounds.append(bound)
        mixed[cost] = stage
    return dict(two_swap_traces=traces, prior_to_cross_bracket=sum(n for _, n in traces), terminal_dimension=sum(low.values()),
                terminal_singlets=trivial_multiplicity(low, "G2"), terminal_weight_ranks=ranks,
                paired_colour_multiplicities=channels, mixed_jet_ranks=mixed, successive_upper_bounds=bounds,
                upper=bounds[-1])


def calculate():
    weight = (6, 6, 6)
    diagrams = [partner(weight, order) for order in (2, 3, 4)]
    assert diagrams == [(5, 5, 3, 3), (6, 6, 3), (7, 7)] and partner(weight, 5) is None
    lower = second_bound()
    assert lower["two_swap_traces"][0][1] == 79 and lower["successive_upper_bounds"] == [100, 39, 9, 9]
    third = projection(6, (3,), 13, True, "G2", 1)[0]
    assert third == 1
    alternative = (5, 5, 5)
    assert kostka(diagrams[1], alternative) == 1 and occupations(alternative, 13) == []
    assert occupations((7, 7), 12) == [(((7, 6), (7, 6)), 0)]
    fourth, traces = projection(7, (), 12, True, "G2")
    assert fourth == 0
    cochains = []
    for order, diagram in zip((2, 3, 4), diagrams):
        archived = json.loads((PRIMARY / "5691072" / f"order{order}.json").read_bytes())
        assert archived["shape"] == list(diagram)
        cochains.append(archived["cochains_by_length"][str(16 - order)])
    return dict(group="G2", source_flavour=weight, source_length=15, diagrams=diagrams,
                degrees_and_colour_dimensions=sorted([e, y, sum(r.values())]
                                                     for (e, y), r in axis("G2").items()),
                target_cochains=cochains, classical_target_upper_bounds=[lower["upper"], 0, fourth],
                order2=lower, order3_single_weight_bound=third,
                order3_alternative_weight=alternative, order4_traces=traces)
