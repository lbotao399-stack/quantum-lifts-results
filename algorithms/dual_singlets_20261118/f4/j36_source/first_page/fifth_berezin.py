"""Integrate the fourteen lifted modes of the glueball flavour at linear order in the fields."""
from fractions import Fraction
from itertools import combinations
import json
from math import comb, factorial, prod

from colour import wedge
from fifth_instanton import calculate


def integral():
    endpoint = calculate()
    colours = endpoint["doublets"]
    pairs = [(colours.index(i), colours.index(j), 2*b) for i,j,b in endpoint["alternating_form"]]
    all_modes = (1 << len(colours)) - 1
    ordered = [i for a,b,_ in pairs for i in (a,b)]
    pfaffian = wedge(ordered)[1] * prod(a for _,_,a in pairs)
    assert pfaffian == -1024

    def power(k):
        rows = {}
        for chosen in combinations(pairs,k):
            key, sign = wedge([i for a,b,_ in chosen for i in (a,b)])
            rows[key] = sign * factorial(k) * prod(a for _,_,a in chosen)
        return rows

    def remove(mask, i):
        if not mask & (1 << i):
            return mask, 0
        return mask ^ (1 << i), (-1)**((mask & ((1 << i)-1)).bit_count())

    def contract(rows):
        result = {}
        for mask, v in rows.items():
            for i,j,a in pairs:
                tail, s = remove(mask,i)
                target, t = remove(tail,j)
                if s*t:
                    result[target] = result.get(target,0) + Fraction(s*t*v,a)
        return result

    old = {wedge([colours.index(i) for i in word])[0]:v for word,v in endpoint["fifth_power"]}
    assert power(5) == old
    checks = []
    for m in (2,3):
        factor = Fraction(factorial(7-m)*pfaffian,factorial(m))
        left = power(7-m)
        hits = 0
        for indices in combinations(range(14),2*m):
            mask, _ = wedge(indices)
            complement = all_modes ^ mask
            dual = left.get(complement,0)
            sign = wedge([i for i in range(14) if complement & (1<<i)] + list(indices))[1]
            contracted = {mask:1}
            for _ in range(m):
                contracted = contract(contracted)
            assert factor * contracted.get(0,0) == sign * dual
            hits += bool(dual)
        assert hits == comb(7,m) and factor.denominator == 1
        checks.append(dict(remainder_degree=2*m, basis_monomials=comb(14,2*m),
                           nonzero_pairings=hits, contraction_coefficient=int(factor)))
    top = power(7)[all_modes]
    assert top == factorial(7)*pfaffian == -5160960
    return dict(doublet_order=colours, pfaffian=pfaffian, top_coefficient=top, tests=checks,
                with_fifth_power_binomial=[checks[0]["contraction_coefficient"],
                                           5*checks[1]["contraction_coefficient"]],
                physical_instanton_supercharge_pairing="open")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(integral(),indent=2))
