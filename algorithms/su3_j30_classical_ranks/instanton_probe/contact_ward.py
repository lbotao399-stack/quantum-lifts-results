"""Apply four exact classical boundaries to all pure four-fermion colours and the mixed test tensor."""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import gcd, lcm
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import runpy

from cartan_contacts import cartan_scalar_expansion
from cubic_colours import source
from four_fermion_contacts import (
    canonical,
    contraction,
    contracts,
    difference,
    independent,
    polarize,
    projected,
    projected_basis,
)
from mixed_contact import contact_image
from source_fluctuations import sign
from verify import HASHES, blob

LAYOUTS = ("000fabcde", "0afbd0ce0", "f0abd/00ce", "0aed/f0b0c")


def field_strength_variation():
    script = blob("dot-r007-t2-primary/replay.py")
    assert sha256(script).hexdigest() == HASHES["replay.py"]
    with TemporaryDirectory() as temp:
        path = Path(temp) / "replay.py"
        path.write_bytes(script)
        kernel = runpy.run_path(str(path))["q_letter"]
        rules = kernel((0, 0, (0, 1, 2)))
    masks = [(0,), (1,), (2,), (1, 2), (0, 2), (0, 1)]
    result = {}
    for coefficient, letters in rules:
        assert all((p, m) == (0, 0) for p, m, _ in letters)
        word = tuple(masks.index(t) for _, _, t in letters)
        result[word] = coefficient * (-1) ** word.count(4)
    assert result == {(i, 3 + i): -1 for i in range(3)} | {
        (3 + i, i): 1 for i in range(3)
    }
    return result


def boundary(layout, variation):
    # 0 is phi^1. The three a,b,c slots are epsilon contracted fermions.
    # d,e contract a matter fermion and scalar; f is the bare triple jet.
    terms = Counter()
    for order in permutations((3, 4, 5)):
        parity = (-1) ** sum(a > b for k, a in enumerate(order) for b in order[k + 1 :])
        for flavour in range(3):
            letters = dict(zip("abcde0f", (*order, 3 + flavour, flavour, 0, 6)))
            traces = tuple(tuple(letters[x] for x in tr) for tr in layout.split("/"))
            j = next(i for i, w in enumerate(traces) if 6 in w)
            w = traces[j]
            k = w.index(6)
            previous = sum(x >= 3 for tr in traces[:j] for x in tr) + sum(
                x >= 3 for x in w[:k]
            )
            for pair, c in variation.items():
                word = w[:k] + pair + w[k + 1 :]
                terms[traces[:j] + (word,) + traces[j + 1 :]] += (
                    parity * (-1) ** previous * c
                )
    return [(c, w) for w, c in terms.items() if c]


def adjoint_input(polynomial, functional):
    # Transpose the prescribed Cartan coordinate onto actual four-derivative
    # colours. A boundary can have three or four off-diagonal differentiated entries.
    totals = Counter()
    weights = {s for s, _ in functional}
    for (s, mask), value in polynomial.items():
        if s not in weights:
            continue
        slots = [i for i in range(24) if mask & (1 << i)]
        for stay in slots:
            if stay % 8 not in (0, 4):
                continue
            cut = [i for i in slots if i != stay]
            counts = Counter(i // 8 for i in cut)
            if sorted(counts.values()) != [1, 1, 2]:
                continue
            doubled = next(i for i, n in counts.items() if n == 2)
            for last in (i for i in cut if i // 8 == doubled):
                ordered = [i for i in cut if i != last] + [last]
                entries = [i % 8 for i in ordered]
                scale = 3 ** sum(i in (0, 4) for i in entries)
                duals = [
                    11 + i // 4 if i in (0, 4) else 3 * (i % 3) + i // 3
                    for i in entries
                ]
                tail = mask
                orientation = 1
                for i in ordered:
                    tail ^= 1 << i
                    orientation *= sign(1 << i, tail)
                for z in range(2):
                    bit = 1 << (8 * doubled + 4 * z)
                    if bit & tail:
                        continue
                    pairing = functional.get((s, tail | bit), 0)
                    if pairing:
                        totals[tuple(duals + [9, 10, 9 + z])] += (
                            value * orientation * sign(bit, tail) * pairing / scale
                        )
    return {point: v for point, v in totals.items() if v}


def constraint(polynomial, mixed, functional, circuits):
    dual = adjoint_input(polynomial, functional)
    values = [
        sum(contraction(circuit, point) * c for point, c in dual.items())
        for circuit in circuits
    ]
    values.append(sum(c * mixed.get(key, 0) for key, c in functional.items()))
    return {j: v for j, v in enumerate(values) if v}


def primitive_relation(rows):
    basis = {}
    for i, row in enumerate(rows):
        v, combo = dict(row), {i: Fraction(1)}
        for pivot, (earlier, weights) in sorted(basis.items()):
            if pivot in v:
                factor = v[pivot]
                difference(v, earlier, factor)
                difference(combo, weights, factor)
        if not v:
            den = lcm(*(c.denominator for c in combo.values()))
            integers = [int(den * combo.get(j, 0)) for j in range(len(rows))]
            common = gcd(*integers)
            return [v // common for v in integers]
        pivot = min(v)
        scale = Fraction(1, v[pivot])
        basis[pivot] = (
            {k: scale * c for k, c in v.items()},
            {k: scale * c for k, c in combo.items()},
        )
    raise ValueError("All Ward rows are independent")


def verify():
    prior, _, target, functional = projected_basis()
    patterns = [p for p in contracts() if projected(p) and p == min(projected(p))]
    circuits = [polarize(projected(p)) for p in patterns]
    labels = [
        patterns.index(canonical(w)) for w in prior["independent_source_patterns"]
    ] + [len(circuits)]
    terms = [(Fraction(c), w) for c, w in source()]
    polynomial, _ = cartan_scalar_expansion(terms)
    mixed = contact_image(terms)
    probes = sorted(functional)
    rows = [constraint(polynomial, mixed, {key: 1}, circuits) for key in probes]
    assert [[str(r.get(i, 0)) for i in labels[:3]] for r in rows] == [
        r[:3] for r in prior["distinguished_coefficients"]
    ]
    rhs = [target.get(key, 0) for key in probes]
    rule = field_strength_variation()
    for layout in LAYOUTS:
        terms = boundary(layout, rule)
        polynomial, _ = cartan_scalar_expansion(terms)
        rows.append(constraint(polynomial, contact_image(terms), functional, circuits))
        rhs.append(0)
    variable_count = len(circuits) + 1
    rank = len(independent(rows))
    augmented = len(
        independent([r | {variable_count: v} if v else r for r, v in zip(rows, rhs)])
    )
    relation = primitive_relation(rows)
    assert all(
        sum(c * r.get(j, 0) for c, r in zip(relation, rows)) == 0
        for j in range(variable_count)
    )
    discrepancy = sum(c * v for c, v in zip(relation, rhs))
    assert (variable_count, rank, augmented) == (87, 7, 8) and discrepancy
    seed_fit = [
        Fraction(c) for c in ("62453/347490", "8071/694980", "-10213/463320", "37/702")
    ]
    residual = sum(c * rows[4].get(i, 0) for i, c in zip(labels, seed_fit))
    assert residual == Fraction(6032431, 45714240)
    return dict(
        boundary_words=list(LAYOUTS),
        tested_colours=variable_count,
        rank=rank,
        augmented_rank=augmented,
        primitive_left_relation=relation,
        target_values=list(map(str, rhs)),
        left_relation_on_target=str(discrepancy),
        first_boundary_on_four_displayed_colours=[
            str(rows[4].get(i, 0)) for i in labels
        ],
        first_boundary_on_fitted_combination=str(residual),
        rows=[{str(i): str(c) for i, c in sorted(r.items())} for r in rows],
    )


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions")
    print(json.dumps(verify(), indent=2))
