"""Solve a full scalar-degree restriction of the SU(3) first equation on GL(2)."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations, product
from math import lcm
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from su3_level32_centralizer import evaluate, restrict_words
from su3_level32_content import incident
from su3_level32_layers import DATA, content, first_targets, formal, joint, physical, residual
from verify import add_scaled, exterior_sign, multiply


def single(slot):
    if formal.ODD[slot // 4]:
        return 0, 1 << slot
    return 1 << (4 * slot), 0


@lru_cache(maxsize=50000)
def change_basis(power, mask):
    """Write A=[[a+h,e],[f,a-h]] in the original order of exterior generators."""
    if not power and not mask:
        return {(0, 0): 1}
    if mask:
        slot = mask.bit_length() - 1
        previous = change_basis(power, mask ^ (1 << slot))
    else:
        slot = (power.bit_length() - 1) // 4
        previous = change_basis(power - (1 << (4 * slot)), mask)
    start, entry = 4 * (slot // 4), slot % 4
    if entry in (0, 3):
        values = {single(start): 1, single(start+1): 1 if entry == 0 else -1}
    else:
        values = {single(start+entry+1): 1}
    return multiply(previous, values)


def separate(row):
    blocks = defaultdict(dict)
    for (power, mask), coefficient in row.items():
        centre = bytearray()
        central_mask = central_power = 0
        for k in range(1, len(formal.ALPH)):
            bit = (mask >> (4*k)) & 1
            exponent = (power >> (16*k)) & 15
            centre.extend([k] * (bit + exponent))
            central_mask |= bit << (4*k)
            central_power |= exponent << (16*k)
        active = power - central_power, mask ^ central_mask
        sign = exterior_sign(central_mask, active[1])
        add_scaled(blocks[bytes(centre)], {active: coefficient}, sign)
    return blocks


def allowed_targets(level):
    if level == 32:
        return first_targets()
    from g2_full_classical_screen_v1 import contents
    return sorted(bytes(sorted(formal.ID[k] for k in t)) for t in contents((1,1,7,4,3), 11))


def first_residual(level, data, vector, keep):
    if level == 32:
        return residual(data, vector, keep)
    rows = {t: {} for t in allowed_targets(level) if keep(t)}
    digest = hashlib.sha256()
    count = 0
    for w, coefficient in joint.s.stream_rows(data/"negative_Q1.jsonl.gz", digest):
        t = content(w)
        if t in rows:
            rows[t][w] = -Fraction(coefficient)
        count += 1
    assert count == 517468 and digest.hexdigest() == formal.INPUT
    for w, coefficient in vector.items():
        for v, n in formal.q(w).items():
            t = content(v)
            if t in rows:
                add_scaled(rows[t], {v: n}, coefficient)
    return rows


def restricted_target(level, data, count):
    targets = first_residual(level, data, {}, lambda t: t.count(0) == count)
    denominator = lcm(*(c.denominator for row in targets.values() for c in row.values()))
    blocks = defaultdict(dict)
    checks = {}
    for t, row in targets.items():
        raw = evaluate(restrict_words(row, denominator))
        split = {}
        for key, value in raw.items():
            add_scaled(split, change_basis(*key), value)
        assert bool(split) == bool(raw)
        checks[t.hex()] = len(raw)
        for centre, active in separate(split).items():
            assert not blocks[centre].keys() & active.keys()
            blocks[centre].update(active)
    return denominator, checks, blocks


@lru_cache(maxsize=None)
def invariants(labels):
    """Spanning products of traceless 2x2 pair traces and at most one triple trace."""
    if not labels:
        return ((),)
    words = set()
    # With pairs only, the first label fixes one pair. The unique triple
    # in odd degree can occur anywhere.
    subsets = (combinations(range(len(labels)), 3) if len(labels) % 2
               else ((0, i) for i in range(1, len(labels))))
    for slots in subsets:
        trace = bytes(labels[i] for i in slots)
        rest = bytes(k for i, k in enumerate(labels) if i not in slots)
        for tail in invariants(rest):
            phase, word = formal.multi((trace,) + tail)
            if phase:
                words.add(word)
    return tuple(sorted(words))


@lru_cache(maxsize=16000)
def active_trace(trace):
    terms = Counter()
    for walk in product(range(2), repeat=len(trace)):
        sign = 1
        power = mask = 0
        for k, i, j in zip(trace, walk, walk[1:]+walk[:1], strict=True):
            slot = 4*k + (1 if i == j else 2+i)
            even, odd = single(slot)
            sign *= exterior_sign(mask, odd)
            if i == j == 1:
                sign = -sign
            if not sign:
                break
            power += even
            mask |= odd
        if sign:
            terms[power, mask] += sign
    return {k: v for k, v in terms.items() if v}


@lru_cache(maxsize=30000)
def active_word(word):
    row = {(0, 0): 1}
    for tr in word:
        row = multiply(row, active_trace(tr))
    return row


def remove_centre(content, centre):
    remainder = list(content)
    for k in centre:
        if k not in remainder:
            return None
        remainder.remove(k)
    return bytes(remainder)


def image(word, centre):
    sign = (-1)**sum(formal.ODD[k] for k in centre)
    row = {}
    for output, coefficient in formal.q(word).items():
        if not any(0 in tr for tr in output):
            add_scaled(row, active_word(output), sign * coefficient)
    return row


def solve_block(centre, target, targets, parents):
    # Different active contents have disjoint monomial coordinates. Choose
    # determining coefficients in each short SU(2) invariant space first.
    selected = []
    for labels in targets:
        active = remove_centre(labels, centre)
        if active is None:
            continue
        rows = [active_word(w) for w in invariants(active)]
        columns = sorted(set().union(*(r.keys() for r in rows)))
        if columns:
            matrix = np.array([[row.get(k, 0) for k in columns] for row in rows], dtype=np.int64)
            selected.extend(columns[i] for i in joint.piv(csr_matrix(matrix)))
    assert len(selected) == len(set(selected))
    words = sorted({w for labels in parents if (active := remove_centre(labels, centre)) is not None
                    for w in invariants(active)})
    images = [image(w, centre) for w in words]
    matrix = csr_matrix(np.array([[r.get(k, 0) for k in selected] for r in images], dtype=np.int64))
    rhs = [-target.get(k, 0) for k in selected]
    if not matrix.nnz:
        if not any(rhs):
            assert not target
            return dict(kind="zero", parents=len(words))
        index = next(i for i, v in enumerate(rhs) if v)
        result, coeffs = dict(kind="exact_separator", image_rank_mod_prime=0), {index: Fraction(1)}
    else:
        result, coeffs = physical.certify(matrix, rhs, lambda _: None)
    result.update(parents=len(words), selected_coordinates=len(selected), rhs_monomials=len(target))
    if result["kind"] == "exact_solution":
        check = target.copy()
        for i, value in coeffs.items():
            add_scaled(check, images[i], value)
        assert not check, "correction fails on omitted polynomial coordinates"
        result["terms"] = [[joint.s.hexes(words[i]), str(v)] for i, v in sorted(coeffs.items())]
    elif result["kind"] == "exact_separator":
        functional = {selected[i]: v for i, v in coeffs.items()}
        assert all(sum(v * row.get(k, 0) for k, v in functional.items()) == 0 for row in images)
        pairing = sum(v * target.get(k, 0) for k, v in functional.items())
        assert pairing
        result["functional"] = [[str(k[0]), str(k[1]), str(v)] for k, v in functional.items()]
        result["pairing_with_scaled_Q1"] = str(pairing)
    else:
        raise ValueError(result)
    return result


def run(level, data, count, out):
    denominator, checks, blocks = restricted_target(level, data, count)
    physical_targets = sorted(t for t in allowed_targets(level) if t.count(0) == count)
    assert physical_targets
    parents = {g for t in physical_targets for g in incident(t) if g.count(0) == count}
    strip = lambda t: bytes(k for k in t if k)
    targets = sorted(map(strip, physical_targets))
    parents = sorted(map(strip, parents))
    records = {}
    for centre, target in sorted(blocks.items(), key=lambda part: (-len(part[0]), part[0])):
        if target:
            records[centre.hex()] = solve_block(centre, target, targets, parents)
            print(centre.hex(), records[centre.hex()]["kind"], len(target), flush=True)
    result = dict(level=level, scalar_degree=count, target_contents=len(targets), parent_contents=len(parents),
                  integer_rhs_denominator=denominator, original_centralizer_monomials=checks,
                  central_coefficients=records)
    joint.s.write(out, result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level", type=int, choices=(32,34), default=32)
    p.add_argument("--data", type=Path)
    p.add_argument("--scalar-degree", type=int, choices=range(2, 8), default=3)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("exact coefficient verification requires enabled assertions")
    data = args.data or (DATA if args.level == 32 else Path(__file__).resolve().parents[1]/"su3_j34_20261118")
    run(args.level, data, args.scalar_degree, args.out)
