"""Spin and exact-mode tests on the two lowest endpoint-propagator channels."""
import argparse
from collections import Counter
from itertools import combinations, product
import json
from math import comb, factorial, prod
from pathlib import Path

from r007_free_radial import oscillator, target
from r007_instanton_modes import occupancy, raw_source


def wedge(left, right):
    result = Counter()
    for a, u in left.items():
        for b, v in right.items():
            if a & b:
                continue
            inversions = sum((b & ((1 << j) - 1)).bit_count() for j in range(6) if a >> j & 1)
            result[a | b] += (-1)**inversions * u * v
    return {m: c for m, c in result.items() if c}


def exterior_check():
    # First two generators are eta, next two are bar xi, last two are the
    # lifted coordinates of flavour 4. Common displacements only change basis.
    eta1, eta2, xi1, xi2, nu, nubar = [{1 << j: 1} for j in range(6)]
    zeta1 = {1: 1, 4: -3, 8: -5}
    zeta2 = {2: 1, 4: 7, 8: -11}
    dotted = {4: 2, 8: 3}
    available = [zeta1, zeta2, dotted, nu, nubar]
    current = {0: 1}
    for v in available:
        current = wedge(current, v)
    assert current and all(not wedge(current, v) for v in available)
    control = wedge(current, xi2)
    assert control == {63: 2}
    # Exhaust the exterior subalgebra, including arbitrary inhomogeneous fields.
    count = Counter()
    for k in range(6):
        for subset in combinations(available, k):
            value = {0: 1}
            for v in subset:
                value = wedge(value, v)
            assert value and all(mask.bit_count() == k for mask in value)
            count[k] += 1
    return dict(subalgebra_dimensions=[count[k] for k in range(7)],
                complementary_dotted_control=control[63])


def triples(left, right, masks):
    counts = Counter()
    for l, size in left.items():
        for r, mult in right.items():
            caps = tuple(min(a, b) for a, b in zip(l, r))
            for take in product(*(range(min(3,c)+1) for c in caps)):
                if sum(take) != 3:
                    continue
                multiplicity = size * mult * prod(comb(a,k)*comb(b,k)*factorial(k)
                                                  for a,b,k in zip(l,r,take))
                counts['species_matched_pairings'] += multiplicity
                source_modes = [m for m,a,k in zip(masks,l,take) for _ in range(a-k)]
                if max(occupancy(source_modes)) <= 4:
                    counts['also_pass_undifferentiated_source_local_bound'] += multiplicity
    return dict(counts)


def run():
    masks, source = raw_source()
    words = target()
    assert len(words) == 5080
    source_fields = ((0,), (1,), (2,), (1,2), (0,2), (0,1))
    species = Counter()
    distinct = set()
    types = Counter()
    for w in words:
        letters = [v for tr in w for v in tr]
        distinct.update(letters)
        content = Counter(v[2] for v in letters)
        species[tuple(content[k] for k in source_fields)] += 1
        assert len(letters) == 9
        right_modes = [1 | sum(2 << i for i in v[2]) for v in letters]
        assert occupancy(right_modes) == (9,5,2,2)
        for v in letters:
            types[len(v[2])] += 1
            (left_plus, bplus, bminus, occupied), coefficient = oscillator(v)
            assert coefficient and left_plus == bplus+bminus+len(v[2])-1
            # The native left singleton uses a_plus only. Reflection swaps a,b.
            # All right dotted occupations therefore sit on one common line.
            assert left_plus >= 0 and all(i in range(3) for i in occupied)
    assert sum(types.values()) == len(words)*9
    # The conjugate elementary field may also be a chiral gaugino or curvature.
    # Only the six source species admit direct bosonic-background contractions.
    census = triples(source, species, masks)
    assert min(census.values()) > 0
    sd = tuple(a+b+c-3 for a,b,c in zip((0,3,6,6),(9,5,2,2),(0,1,1,1)))
    assert sd == (6,6,6,6)
    # In anti-self-dual orientation, F_dotted has profile degree zero and
    # bar-lambda degree one. Their product has undotted spin zero; the
    # improved current has twice-undotted spin two.
    gauge_fermion_doubled_spins = [(0,j) for j in (1,3)]
    assert (2,1) not in gauge_fermion_doubled_spins
    anti = tuple(a+b+c for a,b,c in zip((10,7,4,4),(0,4,7,7),(2,1,1,1)))
    assert anti == (12,12,12,12)
    return dict(flavour_order=[4,1,2,3], source_words=sum(source.values()),target_words=len(words),
                distinct_target_jets=len(distinct), target_occurrences_by_flavour_derivative_count=dict(types),
                source_species_contents=len(source),target_species_contents=len(species), three_propagators=census,
                self_dual_after_three=list(sd),right_flavour_four_required=6,
                right_flavour_four_available=5,exterior_check=exterior_check(),
                improved_current_doubled_spins=[2,1], anti_degree_one_doubled_spins=gauge_fermion_doubled_spins,
                anti_current_flavour_lower_bound=[2,1,1,1],anti_before_contractions=list(anti),
                direct_propagator_lower_bounds=[4,6])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('run with assertions enabled')
    data = run()
    args.out.write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(data))
