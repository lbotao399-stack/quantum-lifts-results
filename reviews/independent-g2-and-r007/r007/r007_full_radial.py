"""Free radial norm of the charge-conjugation-odd SU(3) outgoing candidate.

The trace form and singleton metric are those of r007_free_radial.py.
The calculation treats all nine colour slots before imposing finite-rank identities.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import factorial, prod
from pathlib import Path

import r007_free_radial as radial
import r007_primary_replay as jets


def permutation_sign(order, parity):
    odd = [i for i in order if parity[i]]
    return (-1) ** sum(i > j for k, i in enumerate(odd) for j in odd[k+1:])


def topology(traces):
    pairs = tuple(sorted(tuple(sorted(tr)) for tr in traces if len(tr) == 2))
    cubic, = [tuple(tr) for tr in traces if len(tr) == 3]
    return pairs, tuple(sorted(cubic))


def by_profile(state):
    groups = defaultdict(Counter)
    for traces, coefficient in state.items():
        lengths = tuple(map(len, traces))
        assert sorted(lengths) == [2, 2, 2, 3]
        # Conjugation of every generator supplies (-1)^9. Reversing the
        # graded traces accounts for all signs of the odd field components.
        reversal = (-1)**sum(len(tr) + sum(map(jets.odd,tr))*(sum(map(jets.odd,tr))-1)//2
                               for tr in traces)
        sign, rev = jets.canonical(tuple(tr[::-1] for tr in traces))
        assert state[rev] == -sign * reversal * coefficient
        letters = tuple(x for tr in traces for x in tr)
        inds = sorted(range(9), key=letters.__getitem__)
        slots = {old:new for new,old in enumerate(inds)}
        cut, labels = 0, []
        for n in lengths:
            labels.append(tuple(slots[j] for j in range(cut, cut+n)))
            cut += n
        profile = tuple(letters[i] for i in inds)
        groups[profile][topology(labels)] += coefficient * permutation_sign(inds, tuple(map(jets.odd, letters)))
    return groups


def matchings(profile):
    blocks = [tuple(i for i,x in enumerate(profile) if x==letter) for letter in sorted(set(profile))]
    parities = tuple(map(jets.odd, profile))
    for images in product(*(permutations(b) for b in blocks)):
        perm = tuple(i for group in images for i in group)
        yield perm, permutation_sign(perm, parities)


@lru_cache(None)
def loops(first, second):
    """Count closed adjoint lines; a self-contraction of a cubic gives zero."""
    # Merge each quadratic trace into one adjoint line on the two matched words.
    parent = list(range(9))
    def root(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i,j in first[0] + second[0]:
        parent[root(i)] = root(j)
    ends1 = {root(i) for i in first[1]}
    ends2 = {root(i) for i in second[1]}
    if len(ends1) != 3 or ends1 != ends2:
        return None
    return len({root(i) for i in range(9)} - ends1)


def profile_polynomial(profile, diagrams):
    # Slots of the second word are relabelled by the Wick bijection. For an
    # orthonormal generator basis, each quadratic trace is the identity.
    items = [(d,c) for d,c in diagrams.items() if c]
    acc = Counter()
    for perm, sign in matchings(profile):
        for q, value in items:
            relabelled = topology([tuple(perm[i] for i in tr) for tr in q[0]+(q[1],)])
            for p, other in items:
                n = loops(p, relabelled)
                if n is not None:
                    acc[n] += sign * value * other
    weight = prod(map(radial.inner_letter, profile))
    return {n:k*weight for n,k in acc.items() if k}


def evaluate(polynomial, dimension, cubic_norm):
    return cubic_norm * sum(c*dimension**n for n,c in polynomial.items())


def run():
    state = radial.target()
    groups = by_profile(state)
    total = Counter()
    grams = []
    for profile, terms in sorted(groups.items()):
        coefficient = profile_polynomial(profile, terms)
        grams.append(evaluate(coefficient,8,Fraction(20,3)))
        assert grams[-1] >= 0
        total.update(coefficient)
    # For SU(N), the norm of sym Tr(T_i T_j T_k) is
    # (N^2-1)(N^2-4)/(2N). On its Cartan it is (N-1)(N-2)/N.
    cartan = evaluate(total,2,Fraction(2,3))
    old = json.loads(Path(__file__).with_name('r007-free-radial.json').read_text())
    assert cartan == Fraction(old['Cartan_Gram_defining_trace_unit_singleton'])
    full = sum(grams)
    assert full == evaluate(total,8,Fraction(20,3)) and full > cartan > 0
    return dict(trace_terms=len(state), singleton_profiles=len(groups), nonzero_norm_profiles=sum(x>0 for x in grams),
                topology_coefficients={str(n):str(v) for n,v in sorted(total.items()) if v},
                defining_trace_SU3_symmetric_cubic_norm='20/3', Cartan_symmetric_cubic_norm='2/3',
                full_SU3_free_norm=str(full), Cartan_free_norm=str(cartan))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:
        raise RuntimeError('run with assertions enabled')
    receipt=run()
    args.out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(receipt,flush=True)
