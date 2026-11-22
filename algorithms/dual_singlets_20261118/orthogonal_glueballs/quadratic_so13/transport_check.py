"""Project transported trace identities to the flavour representation of the physical target."""

from collections import Counter
from itertools import product
import gzip
import json
from auxiliary import HERE
from verify import words
from terminal_commutators import parity
from functools import lru_cache
from fractions import Fraction
from math import factorial, prod
from flint import fmpq_mat
from second import substitution_on_section
from orthogonal_glueball_relation_span import determinant
from joint import base_relation, affected_sites


@lru_cache(None)
def E(w, i, j):
    out = Counter()
    for pos, t in enumerate(w):
        if not t[2 + j]:
            continue
        v = list(t)
        v[2 + i] += 1
        v[2 + j] -= 1
        raw = w[:pos] + (tuple(v),) + w[pos + 1 :]
        odd = [key(x) for x in raw if parity(x)]
        if len(set(odd)) == len(odd):
            sg = (-1) ** sum(a > b for k, a in enumerate(odd) for b in odd[k + 1 :])
            out[tuple(sorted(raw, key=key))] += sg * t[2 + j]
    return {v: c for v, c in out.items() if c}


def act(poly, i, j):
    new = Counter()
    for w, c in poly.items():
        for v, n in E(w, i, j).items():
            new[v] += n * c
    return {v: c for v, c in new.items() if c}


def project(poly):
    # Use the gl(3) Casimir on the fixed weight (4,2,2).
    mu = (4, 2, 2)

    def quadratic(w):
        return sum(x * (x + 2 - 2 * i) for i, x in enumerate(w))

    vals = set()
    for a in range(4, 9):
        for b in range(a + 1):
            c = 8 - a - b
            t = (a, b, c)
            if b >= c >= 0 and a + b >= 6 and t != mu:
                val = quadratic(t) - quadratic(mu)
                assert val > 0 and val % 2 == 0
                vals.add(val // 2)
    for val in sorted(vals):
        result = Counter({w: c * val for w, c in poly.items()})
        for i, j in ((0, 1), (1, 2), (0, 2)):
            result.subtract(act(act(poly, i, j), j, i))
        poly = {w: c for w, c in result.items() if c}
    assert not act(poly, 0, 1) and not act(poly, 1, 2)
    return {w: Fraction(c, prod(vals)) for w, c in poly.items()}


def transported():
    base = HERE.parent / "scalar_multipliers"
    basis = words((5, 5, 4, 2, 2))
    old = words((4, 4, 4, 2, 2))
    # Reconstruct the first two families from the public ingredients.
    from second import pure_relation, pure_specs

    r7 = json.loads(
        gzip.decompress((base.parent / "so11_j36" / "relations13.json.gz").read_bytes())
    )
    rows = []
    for record in r7["independent"]:
        rel = determinant(7, record["spec"])
        pol = {
            w: c for w in words((5, 5, 2, 2, 2)) if (c := rel.commuting_coefficient(w))
        }
        for letter in range(5):
            rows.append(substitution_on_section(pol, letter=letter))
        rel.clear()
    specs0 = list(pure_specs())
    for j in (0, 2, 28):
        rel = pure_relation(specs0[j])
        pol = {
            w: c for w in words((6, 6, 2, 2, 2)) if (c := rel.commuting_coefficient(w))
        }
        rows.append(
            substitution_on_section(substitution_on_section(pol, spin=1), spin=0)
        )
        rel.clear()
    mixed = json.loads((base / "mixed.json").read_text())["rows"]
    for j in (0, 2, 3, 38, 39):
        rel = base_relation(mixed[j])
        p = {w: c for w in old if (c := rel.commuting_coefficient(w))}
        for f in range(3):
            v = Counter()
            for w, spec in affected_sites(old, basis, f):
                for col, m in spec.items():
                    v[basis[col]] += p.get(w, 0) * m
            rows.append({w: c for w, c in v.items() if c})
        rel.clear()
    return rows


def key(v):
    return sum(v), v


@lru_cache(None)
def trace_coords(site_tuple):
    # Restrictions of individual sorted trace monomials to the labelled planes in site_tuple.
    out = Counter()
    for assignments in product(*(words(v) for v in site_tuple)):
        seq = tuple(t for s in assignments for t in s)
        odds = [key(t) for t in seq if parity(t)]
        if len(set(odds)) < len(odds):
            continue
        exp = sum(a > b for i, a in enumerate(odds) for b in odds[i + 1 :])
        for plane in assignments:
            odd = [j for t in plane for j in range(2) if t[j]]
            exp += sum(a > b for i, a in enumerate(odd) for b in odd[i + 1 :])
        n = prod(factorial(k) for k in Counter(seq).values()) // prod(
            factorial(k) for s in assignments for k in Counter(s).values()
        )
        out[tuple(sorted(seq, key=key))] += (
            n * 2 ** len(seq) * (-1) ** (exp + sum(map(sum, site_tuple)) // 2)
        )
    return {k: v for k, v in out.items() if v}


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with assertions enabled")
    rows = transported()
    basis = words((5, 5, 4, 2, 2))
    unprojected = fmpq_mat([[str(row.get(w, 0)) for w in basis] for row in rows]).rank()
    projected = [project(v) for v in rows]
    archive = json.loads(gzip.decompress((HERE / "rows.json.gz").read_bytes()))
    sites = [(i, w) for i, w in enumerate(basis) if len(w) > 6]
    direct = [[r[3].get(str(i), "0") for i, w in sites] for r in archive]
    old = [
        [
            str(sum(row.get(v, 0) * c for v, c in trace_coords(w).items()))
            for i, w in sites
        ]
        for row in projected
    ]
    ranks = [fmpq_mat(x).rank() for x in (direct, old, direct + old)]
    assert unprojected == 23 and ranks == [20, 15, 20]
    print("Transported rank before flavour projection", unprojected)
    print(
        "Direct, flavour-projected transported, combined ranks on 953 coordinates",
        ranks,
    )
