"""Resolve rank-two odd BMN source weights by classical ranks and the primary index."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

from generator_bounds import dressings, first_failure, table
from verify import HERE, cochains, read, receipt


G2_JOBS = ("5691136", "5691142", "5691420", "5691422", "5691423", "5691483", "5691485", "5691572")
USP4_JOBS = ("7658772", "7658842", "7658843", "7658854", "7658855", "7658856", "5700336", "5700337")


def columns(group, jobs, maximum):
    _, dimensions = cochains(maximum, group)
    evidence = {}
    for job in jobs:
        receipt(job)
        for name in read(HERE/job/"receipt.json")["files"]:
            row = read(HERE/job/name)
            weight = tuple(row["weight"])
            assert row.get("group", "G2") == group and row["prime"] == 2147483647
            for block in row["blocks"]:
                length = block["length"]
                if length % 2 == 0:
                    continue
                entry = evidence.setdefault(weight, {})
                exact_dimension = dimensions[length][weight][0]
                incoming, outgoing = block["incoming"], block["outgoing"]
                assert 0 <= incoming <= dimensions.get(length-1, {}).get(weight, [0])[0]
                assert 0 <= outgoing <= dimensions.get(length+1, {}).get(weight, [0])[0]
                assert incoming + outgoing <= exact_dimension == block["dimension"]
                assert length not in entry
                entry[length] = dict(job=job, dimension=exact_dimension, incoming=incoming,
                                     outgoing=outgoing, rational_upper=exact_dimension-incoming-outgoing)
    result = []
    for weight, entries in sorted(evidence.items()):
        odd = [y for y, weights in dimensions.items() if y % 2 and weight in weights]
        possible = [y for y in odd if y not in entries or entries[y]["rational_upper"] > 0]
        result.append(dict(weight=weight, ranks=dict(sorted(entries.items())), possible_lengths=possible))
    return result


def symplectic_fourth(columns):
    """Find all third generators that could meet both deficits by level 48."""
    indices = HERE / "rank_two_indices"
    ordinary = table(read(indices / "C2_ordinary_exact.json"))
    primary = table(read(indices / "C2_primary_exact.json"))
    small = json.loads(gzip.decompress((HERE.parent / "b2_level36/odd_through_j28.json.gz").read_bytes()))
    ranks = small["nonempty_blocks"]
    assert len(ranks) == 297
    assert all(y % 2 and sum(w) <= 14 and d == i + o for w, y, d, i, o in ranks)
    targets = ((12,5,5), (8,8,8))
    assert [primary[w] for w in targets] == [-3, -2]

    def capacity(target, weight):
        terms = dressings(target, weight)
        assert all(sum(d) <= 9 for d in terms)
        return sum(n * ordinary[d] for d,n in terms.items())

    old = (10,4,4)
    previous = [2 * capacity(t, old) for t in targets]
    needed = [-primary[t]-n for t,n in zip(targets, previous)]
    assert previous == [2,0] and needed == [1,2]
    possible = []
    for a in range(9):
        for b in range(a+1):
            for c in range(b+1):
                weight = (a,b,c)
                if sum(weight) < 15:
                    continue
                values = [capacity(t, weight) for t in targets]
                if all(x >= n for x,n in zip(values, needed)):
                    possible.append(dict(weight=weight, capacity=values))
    excluded = {r["weight"] for r in columns if not r["possible_lengths"]}
    assert len(possible) == 6 and {r["weight"] for r in possible} == excluded
    # This weight passes every earlier tabulated negative-index inequality.
    last = first_failure({old:2,(8,5,3):1}, ordinary, primary)
    assert last == dict(target=targets[1], level=48, required=2, upper_by_core=[0,1])
    return dict(targets=targets, primary_index=[-3,-2], capacity_of_known_pair=previous,
                required_from_third_generator=needed, excluded_candidates=possible,
                first_forcing_level_in_this_test=48, earlier_index_control=last,
                earlier_index_control_weight=(8,5,3), odd_generator_lower_bound=4)


def compute():
    g2 = columns("G2", G2_JOBS, 18)
    symplectic = columns("USp(4)", USP4_JOBS, 16)
    open_g2 = {r["weight"]: r["possible_lengths"] for r in g2 if sum(r["weight"]) >= 15}
    assert open_g2 == {(6,5,4): [], (6,6,4): [], (6,6,6): [11,13,15]}
    assert all(not r["possible_lengths"] for r in symplectic)
    receipt("7658750")
    ordinary = table(read(HERE/"7658750/G2_ordinary.json"))
    primary = table(read(HERE/"7658750/G2_primary.json"))
    earlier = read(HERE/"odd_through_28.json")["rows"]
    assert len(earlier) == 293
    assert all(length % 2 and sum(weight) <= 14 and dimension == incoming+outgoing
               for weight, length, dimension, incoming, outgoing in earlier)
    first = (6,6,6)
    level44 = (8,7,7)
    reaching_first = [weight for weight in ordinary if sum(weight) >= 15 and
                      sum(n*ordinary[w] for w,n in dressings(first,weight).items())]
    assert sorted(reaching_first) == sorted(open_g2)
    assert (primary[first], primary[level44], primary[19,6,6]) == (-1,-2,-1)
    assert dressings(level44, first) == {(2,1,1): 1} and ordinary[2,1,1] == 1
    assert dressings((19,6,6), first) == {(13,0,0): 1} and ordinary[13,0,0] == 0
    assert first_failure({first: 1}, ordinary, primary) == dict(target=level44, level=44, required=2, upper_by_core=[1])
    candidates = []
    for weight in sorted(ordinary):
        if sum(weight) < 15 or any(w > t for w, t in zip(weight, level44)):
            continue
        if weight in open_g2 and not open_g2[weight]:
            continue
        multiplicities = dressings(level44, weight)
        capacity = sum(n*ordinary[w] for w, n in multiplicities.items())
        if capacity:
            proposal = Counter([first, weight])
            failure = first_failure(proposal, ordinary, primary)
            candidates.append(dict(weight=weight, level=2*sum(weight), level44_capacity=capacity, first_failure=failure))
    remaining = [r["weight"] for r in candidates if r["first_failure"] is None]
    assert len(candidates) == 33 and len(remaining) == 22
    previous = read(HERE/"C2_level192_bounds.json")
    excluded = {r["weight"] for r in symplectic}
    possible_usp = [t["core"] for t in previous["trials"] if t["first_failure"] is None and tuple(t["core"]) not in excluded]
    fourth = symplectic_fourth(symplectic)
    assert not possible_usp
    return dict(G2_columns=g2, USp4_columns=symplectic, G2_required_level36_weight=first,
                G2_odd_generator_lower_bound=2, G2_odd_generators_minimum_distinct_weights=2,
                G2_two_generator_hypothesis=candidates, G2_two_generator_alternatives=remaining,
                USp4_single_additional_generator_alternatives=possible_usp,
                USp4_fourth_generator=fourth)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Run with certificate checks enabled")
    args.out.write_text(json.dumps(compute(), indent=2) + "\n")
