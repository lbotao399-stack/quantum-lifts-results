"""Rebuild the full derivative-free Levi restriction of SU(3) R-007 and its rational primitives."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from io import BytesIO
from itertools import permutations, product
import json
from pathlib import Path
from zipfile import ZipFile

from collision_support import checked_archive_member
from verify import blob

HERE = Path(__file__).resolve().parent

PERMS = tuple(permutations(range(3)))


def phase(xs):
    return (-1) ** sum(u > v for i, u in enumerate(xs) for v in xs[i + 1 :])


def ordered(w):
    ff = [x for x in w if x >= 12]
    if len(set(ff)) != len(ff):
        return (), 0
    return tuple(sorted(w)), phase(ff)


@cache
def orbit(w):
    counts = [sum(x % 3 == i for x in w) for i in range(3)]
    if len({n % 2 for n in counts}) != 1:
        return (), 0
    values = {}
    for perm in PERMS:
        tr = [(x // 3) * 3 + perm[x % 3] for x in w]
        key, sg = ordered(tr)
        sg *= phase(perm) ** (counts[0] % 2)
        if key in values and values[key] != sg:
            return (), 0
        values[key] = sg
    best = min(values)
    return best, values[best]


@cache
def variation(w):
    out = Counter()
    prefix = 1
    for n, x in enumerate(w):
        i, a = divmod(x, 3)
        b, c = (a + 1) % 3, (a + 2) % 3
        terms = []
        if i == 3:
            for f in range(3):
                terms.extend(
                    [(1, 3 * f + b, 12 + 3 * f + c), (-1, 3 * f + c, 12 + 3 * f + b)]
                )
        if i >= 4:
            j, k = (i - 4 + 1) % 3, (i - 4 + 2) % 3
            terms = [(1, 3 * j + b, 3 * k + c), (-1, 3 * j + c, 3 * k + b)]
        for sign, A, B in terms:
            ww, s = ordered(w[:n] + (A, B) + w[n + 1 :])
            if s:
                out[ww] += s * sign * prefix
        if i >= 4:
            prefix = -prefix
    return {k: v for k, v in out.items() if v}


def coinvariant(row):
    out = Counter()
    for w, val in row.items():
        key, sg = orbit(w)
        if sg:
            out[key] += val * sg
    return {k: v for k, v in out.items() if v}


def monomial(e, mask):
    return tuple(i for i, n in enumerate(e) for _ in range(n)) + tuple(
        12 + j for j in range(9) if mask >> j & 1
    )


def stars(n):
    return ((i, j, n - i - j) for i in range(n + 1) for j in range(n - i + 1))


def subtract(x, y, a):
    for k, v in y.items():
        val = x.get(k, 0) - a * v
        if val:
            x[k] = val
        else:
            x.pop(k, None)


def reduce(row, pivots, store):
    expression = {}
    while row:
        key = min(row)
        if key not in pivots:
            break
        i, x, expr = pivots[key]
        factor = row[key]
        subtract(row, x, factor)
        if store:
            subtract(expression, expr, -factor)
    return row, expression


@cache
def tr(word):
    # 0..8 are the SU(2) scalars, 12..20 its fermions; 21..23 are central fermions.
    # The first scalar alone is shifted by h; the central fourth component f is absent in F.
    terms = Counter()
    for options in product(*(range(4) if x == 0 or x >= 3 else range(3) for x in word)):
        last = -1
        iphase = 0
        h = 0
        slots = []
        for x, c in zip(word, options):
            if c == 3:
                if x == 0:
                    h += 1
                else:
                    slots.append(21 + x - 3)
                continue
            slots.append(3 * x + c if x < 3 else 12 + 3 * (x - 3) + c)
            if last == -1:
                last = c
            elif c == last:
                last = -1
            else:
                iphase += 1 if c == (last + 1) % 3 else -1
                last = 3 - last - c
        if last != -1:
            continue
        key, sg = ordered(slots)
        if not sg:
            continue
        trace = 2
        if all(v == 3 for v in options):
            trace += (-2) ** len(options)
        terms[h, key, iphase % 2] += trace * sg * (-1) ** ((iphase % 4) // 2)
    return {w: c for w, c in terms.items() if c}


def multiply(a, b):
    out = Counter()
    for (h, w, i), x in a.items():
        for (k, v, j), y in b.items():
            word, sg = ordered(w + v)
            if sg:
                out[h + k, word, (i + j) % 2] += sg * x * y * (-1) ** (i * j)
    return {w: c for w, c in out.items() if c}


def expand(raw):
    # Group equal products of leading trace factors across the 6,139 terms.
    @cache
    def initial(m):
        if len(m) == 1:
            return tr(m[0])
        return multiply(initial(m[:-1]), tr(m[-1]))

    grouped = {}
    for coefficient, mon in raw:
        f = Fraction(coefficient)
        mon = tuple(tuple(w) for w in mon)
        if len(mon) == 1:
            key = ()
            last = mon[0]
        else:
            key = mon[:-1]
            last = mon[-1]
        bucket = grouped.setdefault(key, Counter())
        for w, c in tr(last).items():
            bucket[w] += f * c
    fit = Counter()
    for prev, row in grouped.items():
        if not prev:
            fit.update(row)
            continue
        fit.update(multiply(initial(prev), {w: c for w, c in row.items() if c}))
    blocks = {}
    for (h, w, iph), c in fit.items():
        if c:
            central = tuple(x - 21 for x in w if x >= 21)
            noncentral = tuple(x for x in w if x < 21)
            b = blocks.setdefault((h, central, iph), {})
            b[noncentral] = c
    return blocks


@cache
def reynolds(w):
    cnt = [sum(x % 3 == c for x in w) for c in range(3)]
    assert len({n % 2 for n in cnt}) == 1
    row = Counter()
    for perm in permutations(range(3)):
        key, s = ordered(tuple(3 * (x // 3) + perm[x % 3] for x in w))
        row[key] += Fraction(s * phase(perm) ** (cnt[0] % 2), 6)
    return {k: v for k, v in row.items() if v}


def expand_reynolds(terms):
    out = Counter()
    for w, c in terms.items():
        for k, v in reynolds(w).items():
            out[k] += c * v
    return {k: v for k, v in out.items() if v}


def act(row):
    out = Counter()
    for w, c in row.items():
        for k, v in variation(w).items():
            out[k] += c * v
    return {k: v for k, v in out.items() if v}


@cache
def pivots(weight, length):
    reps = set()
    raw_count = 0
    for mask in range(512):
        psi = [((mask >> (3 * i)) & 7).bit_count() for i in range(3)]
        for nf in range(min(weight) + 1):
            scal = [k - nf - sum(psi) + psi[i] for i, k in enumerate(weight)]
            if min(scal) < 0 or sum(scal) + sum(psi) + nf != length:
                continue
            for deg in product(*(stars(k) for k in scal + [nf])):
                w = monomial(sum(deg, ()), mask)
                ww, sg = orbit(w)
                raw_count += 1
                if sg:
                    reps.add(ww)
    reps = sorted(reps)
    raw = [coinvariant(variation(v)) for v in reps]
    piv = {}
    for n, v in sorted(enumerate(raw), key=lambda t: len(t[1])):
        row, old = reduce({k: Fraction(c) for k, c in v.items()}, piv, True)
        if row:
            key = min(row)
            q = row[key]
            expr = {j: -c / q for j, c in old.items()}
            expr[n] = 1 / q
            piv[key] = (n, {k: c / q for k, c in row.items()}, expr)
    return reps, piv, raw_count


def test(target, weight):
    length = len(next(iter(target)))
    reps, piv, raw = pivots(tuple(weight), length - 1)
    assert not act(target)
    residual, coeffs = reduce(coinvariant(target), piv, True)
    lifted = {reps[j]: v for j, v in coeffs.items()}
    exact = act(expand_reynolds(lifted))
    difference = {
        w: target.get(w, 0) - exact.get(w, 0)
        for w in target.keys() | exact.keys()
        if target.get(w, 0) != exact.get(w, 0)
    }
    assert coinvariant(difference) == residual
    assert bool(difference) == bool(residual)
    return dict(
        weight=weight,
        incoming_length=length - 1,
        parents=raw,
        octahedral_columns=len(reps),
        octahedral_rank=len(piv),
        target_terms=len(target),
        residual_terms=len(difference),
        primitive=[[list(w), str(v)] for w, v in sorted(lifted.items())],
    )


def rotation(row, axis):
    out = Counter()
    for w, value in row.items():
        for n, x in enumerate(w):
            species, c = divmod(x, 3)
            if c == axis:
                continue
            z = 3 * species + 3 - axis - c
            key, sign = ordered(w[:n] + (z,) + w[n + 1 :])
            if sign:
                out[key] += value * sign * (1 if (c - axis) % 3 == 1 else -1)
    return {k: v for k, v in out.items() if v}


def flavour(row, i, j):
    out = Counter()
    for w, c in row.items():
        for n, x in enumerate(w):
            group, col = divmod(x, 3)
            if group not in (j, 4 + i):
                continue
            new, ref = (3 * i + col, 1) if group == j else (12 + 3 * j + col, -1)
            key, sg = ordered(w[:n] + (new,) + w[n + 1 :])
            if sg:
                out[key] += c * sg * ref
    return {k: v for k, v in out.items() if v}


def inputs():
    with ZipFile(
        BytesIO(blob("dot-r007-exact-input/phase34-public-q0-q1-replay.zip"))
    ) as a:
        source = checked_archive_member(
            a,
            "dmk_o330_hw.json",
            "8bee15e4f741312ca2c05e711b6aec8afcd69fe47dddf78e7671a05f5bab9f34",
        )
    with ZipFile(
        BytesIO(blob("dot-j24-formal-completion/j24-phase7-public-full-chain.zip"))
    ) as a:
        frozen = checked_archive_member(
            a,
            "o0_classical_certificate.json",
            "78c4eab761234f9c136a42c1481a3c92cc03887e8345545c0e6810995e019eda",
        )
    j = json.loads(frozen)
    return json.loads(source)["terms"], {
        monomial(*j["target_basis"][i]): Fraction(c) for i, c in j["O0"]
    }


def verify(rebuild=False):
    raw, o = inputs()
    blocks = expand(raw)
    assert (
        len(raw) == 6139
        and len(o) == 117
        and all(not flavour(o, i, j) for i, j in permutations(range(3), 2))
    )
    frozen = json.loads((HERE / "levi_primitives.json").read_bytes())
    assert len(frozen["central_coefficients"]) == len(blocks) == 16

    def check(target, record):
        assert len(target) == record["target_terms"] and all(
            not rotation(target, i) for i in range(3)
        )
        primitive = {tuple(w): Fraction(c) for w, c in record["primitive"]}
        assert act(expand_reynolds(primitive)) == target
        if rebuild:
            assert test(target, record["weight"]) == record
        return len(primitive)

    terms = 0
    npr = 0
    seen = set()
    for row in frozen["central_coefficients"]:
        key = (row["scalar_power"], tuple(row["central_fermions"]), 1)
        assert key not in seen
        seen.add(key)
        target = blocks[key]
        terms += len(target)
        if key == (3, (), 1):
            assert target == {w: -Fraction(5184, 5) * c for w, c in o.items()}
        else:
            algebra = {
                k: v
                for k, v in row.items()
                if k not in ("scalar_power", "central_fermions")
            }
            npr += check(target, algebra)
    assert seen == set(blocks)
    dressing = Counter()
    for w, c in o.items():
        for colour in range(3):
            dressing[tuple(sorted((*w, colour, colour)))] += c
    ng = check(dressing, frozen["scalar_square_O0"])
    return dict(
        source_words=len(raw),
        nonzero_Levi_coefficients=16,
        raw_monomials=terms,
        exact_coefficients=15,
        primitive_orbits=npr,
        cubic_coefficient_real=0,
        cubic_coefficient_imaginary="-5184/5",
        frozen_operator_terms=len(o),
        square_dressing_terms=len(dressing),
        square_primitive_orbits=ng,
        flavour_shear_residuals=[0] * 6,
        nonperturbative_current_matrix_element="open",
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--rebuild",
        action="store_true",
        help="Also repeat the rational elimination for every primitive.",
    )
    assert __debug__
    print(json.dumps(verify(ap.parse_args().rebuild), indent=2))
