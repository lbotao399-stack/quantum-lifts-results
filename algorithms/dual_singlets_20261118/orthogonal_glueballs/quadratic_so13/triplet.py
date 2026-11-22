"""Antisymmetric scalar projector for the highest dotted-spin triplet."""

from collections import Counter
from functools import lru_cache
import json

from auxiliary import HERE
from verify import determinant, encode, first_duals, path_types, unit, words
from orthogonal_glueball_first_pair import PfaffianCoefficient


@lru_cache(None)
def wedge(order, a, b):
    proj = PfaffianCoefficient(order, 0)
    out = Counter()
    for x, c in proj.weights.items():
        for pair, sg in [((a, b), 1), ((b, a), -1)]:
            weight = x
            n = c
            for i, k in enumerate(pair):
                off = 3 * (3 * i + k)
                n *= 1 + ((weight >> off) & 7)
                weight += 1 << off
            out[weight] += sg * n
    proj.weights = {x: n for x, n in out.items() if n}
    div = set()
    for x in proj.weights:
        s = {0}
        for j in range(9):
            for _ in range((x >> (3 * j)) & 7):
                s |= {v + (1 << (3 * j)) for v in s}
        div.update(s)
    proj.divisors = div
    return proj


def compute(colours, power, flavours):
    target = (power + flavours, power) + tuple(int(i < flavours) for i in range(3))
    assert len(words(target)) == 1
    length = 2 * (power + flavours)
    assert not path_types(power + flavours, power, flavours, length - colours - 1)
    doubles = power + flavours - 4
    parts = ["mp"] * doubles
    spec = dict(
        paths=parts,
        single=[(0, 1, 2), tuple(range(3 - flavours))],
        matrices=[unit(4, 3)] * doubles,
        scalar_directions=[[]] * doubles,
    )
    relation = determinant(power, spec)
    relation.project = wedge(1, 1, 2) if flavours == 2 else PfaffianCoefficient(2, 0)
    relation.letters, relation.pf_degree = length + 1, power + flavours + 5
    assert relation.pf_degree == relation.letters - doubles > colours
    rows, counts = first_duals(colours, target, [0])
    expected = 69 if flavours == 2 else (48 if power == 1 else 81)
    assert counts == [expected]
    result = sum(c * relation.coefficient(encode(w)) for w, c in rows[0].items())
    relation.clear()
    return [colours, power, str(result)]


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    for m, last, name in [(2, 13, "triplet"), (3, 16, "quartet")]:
        results = [
            compute(n, p, m)
            for n in range(7, last + 1)
            for p in (n - 4 - m, n - 3 - m)
            if p + m >= 4
        ]
        assert results == json.loads((HERE / (name + ".json")).read_bytes())
        print(name, json.dumps(results))
