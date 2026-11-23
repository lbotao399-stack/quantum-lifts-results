"""Construct spin singlets at the four-derivative-pair endpoint of SU(3) R-007."""
from collections import Counter
from functools import cache
from itertools import product
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "dual_singlets_20261118"))
from flint import nmod_mat
from mixed33_exact import apply, exact
from rank_two_bmn_product_replay import witness
from sp3_omega_product_replay import odd_jet, odd_slot, ordered
from sp3_omega_products import all_products, factor_polynomial, raising
from su3_bmn_parity_products import labelled


@cache
def diagonal_factor(spec):
    """Evaluate a single trace of any allowed length on diag(x,y,-x-y)."""
    result = Counter()
    length = sum(spec[:5])
    for word, coefficient in factor_polynomial(spec).items():
        result[word, ()] += coefficient
        result[(), word] += coefficient
        for allocation in product((0,1), repeat=length):
            left = tuple(x for x, s in zip(word, allocation) if s == 0)
            right = tuple(x for x, s in zip(word, allocation) if s == 1)
            swaps = sum(odd_jet(x) * odd_jet(y) for i, x in enumerate(word) if allocation[i]
                        for j, y in enumerate(word[i+1:], i+1) if not allocation[j])
            result[left, right] += coefficient * (-1)**(swaps + length)
    return {k: v for k, v in result.items() if v}


@cache
def diagonal_product(recipe):
    if not recipe:
        return {((), ()): 1}
    result = Counter()
    for (left, right), c in diagonal_product(recipe[:-1]).items():
        for (u, v), d in diagonal_factor(recipe[-1]).items():
            s, first = ordered(left + u, odd_jet)
            t, second = ordered(right + v, odd_jet)
            result[first, second] += c*d*s*t*(-1)**(odd_slot(right)*odd_slot(u))
    return {k: v for k, v in result.items() if v}


def full_product_weight(plus, minus):
    families = all_products(7, (plus, minus, 3, 0, 0), 7, include_odd_powers=True)
    recipes = sorted({r for group in families.values() for r in group},
                     key=lambda r: (max(sum(g[:5]) for g in r), r))
    columns = [diagonal_product(r) for r in recipes]
    prefix = sum(max(sum(g[:5]) for g in r) <= 3 for r in recipes)
    assert all(row == labelled(r, True) for r, row in zip(recipes[:prefix], columns))
    independent, _, _ = exact(columns)
    assert max(independent) < prefix
    return dict(recipes=len(recipes), quadratic_and_cubic=prefix, rational_rank=len(independent))


@cache
def spin_raise(monomial):
    """Act directly on divided jets after taking the two-colour Cartan trace."""
    result = Counter()
    for i, slot in enumerate(monomial):
        for j, (plus, minus, flavour) in enumerate(slot):
            if minus:
                changed = (plus + 1, minus - 1, flavour)
                sign, word = ordered(slot[:j] + (changed,) + slot[j + 1:], odd_jet)
                if sign:
                    result[monomial[:i] + (word,) + monomial[i + 1:]] += sign * (plus + 1)
    return tuple((k, v) for k, v in result.items() if v)


def verify():
    families = all_products(3, (4, 4, 3, 0, 0), 7, include_odd_powers=True)
    source, rows = [], []
    for family in sorted(families):
        recipes = families[family]
        raising_images = [raising(recipe) for recipe in recipes]
        for combination in exact(raising_images)[2]:
            assert not apply(combination, raising_images)
            row = Counter()
            for index, coefficient in combination.items():
                for monomial, value in labelled(recipes[index], True).items():
                    row[monomial] += coefficient * value
            row = {k: v for k, v in row.items() if v}
            if not row:
                continue
            raised = Counter()
            for monomial, coefficient in row.items():
                for changed, value in spin_raise(monomial):
                    raised[changed] += value * coefficient
            assert not any(raised.values())
            rows.append(row)
            source.append([[recipes[i], c] for i, c in sorted(combination.items())])
    columns = sorted({k for row in rows for k in row})
    # The field chooses independent rows; the archived determinant is over the integers.
    selection = nmod_mat([[int(r.get(k, 0)) for r in rows] for k in columns], 1000003)
    reduced, rank = selection.rref()
    pivots = [next(j for j in range(len(rows)) if reduced[i, j]) for i in range(rank)]
    integral = witness([rows[i] for i in pivots])
    full = [full_product_weight(4, 4), full_product_weight(5, 3)]
    assert full[0]["rational_rank"] - full[1]["rational_rank"] == integral["rank"] == 33
    return dict(weight=[3, 0, 0, 4, 4], length=7, full_trace_weights=full,
                product_recipes=sum(map(len, families.values())), highest_weight_recipes=[source[i] for i in pivots],
                integer_minor=integral)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), separators=(",", ":")))
