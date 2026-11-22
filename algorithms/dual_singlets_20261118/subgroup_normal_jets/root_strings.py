"""Independent root-string calculation of the two ordered colour fans."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path


def g2():
    # Short-root-first Dynkin coordinates. These are the twelve roots obtained
    # from (2,-1), (-3,2) and the two simple reflections.
    positive = [(2,-1), (-3,2), (-1,1), (1,0), (3,-1), (0,1)]
    roots = {tuple(s*x for x in root) for root in positive for s in (-1,1)}
    grade = lambda root: root[0] % 3
    return roots, grade, lambda r: r[1], lambda r: Fraction(2,3) if grade(r) else 2, 3


def f4():
    # Twice the Cartesian root coordinates. The coroot used for the two
    # external legs is e1-e2, of squared length two in the inherited form.
    roots = set(product((-1,1), repeat=4))
    for n in (1,2):
        for support in combinations(range(4), n):
            for signs in product((-2,2), repeat=n):
                v = [0]*4
                for i,s in zip(support,signs):
                    v[i] = s
                roots.add(tuple(v))
    assert len(roots) == 48
    return (roots, lambda r: r[0]%2, lambda r: Fraction(r[0]-r[1],2),
            lambda r: Fraction(sum(i*i for i in r),4), 2)


def terms(roots, grade, external, square, period):
    first, second = Counter(), Counter()
    zero = (0,)*len(next(iter(roots)))
    for alpha in roots:
        a = grade(alpha)
        first[tuple((a,-a % period).count(k) for k in range(1,period))] += Fraction(external(alpha)**2,2)
        for beta in roots:
            def reach(sign):
                length = 0
                while True:
                    weight = tuple(y + sign*(length+1)*x for x,y in zip(alpha,beta))
                    if weight != zero and weight not in roots:
                        return length
                    length += 1

            below, above = reach(-1), reach(1)
            b = grade(beta)
            slots = a, (b-a)%period, -b%period
            powers = tuple(slots.count(k) for k in range(1,period))
            second[powers] += Fraction(square(alpha)*below*(above+1)*external(alpha)*external(beta),4)
    return [{key:q for key,q in result.items() if q} for result in (first,second)]


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    saved = json.loads(Path(__file__).with_name("colour_fans.json").read_bytes())
    for inputs, row in zip((g2(), f4()), saved, strict=True):
        first, second = terms(*inputs)
        for result, key in ((first,"one_loop_colour"),(second,"two_loop_ordered_colour")):
            assert result == {tuple(item["normal_occupations"]):Fraction(item["coefficient"])
                              for item in row[key]}
        print(row["embedding"], sum(first.values()), sum(second.values()))
