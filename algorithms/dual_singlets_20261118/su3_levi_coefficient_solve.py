"""Test a four-spectator coefficient against all two-colour classical boundaries."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from g2_full_classical_screen_v1 import contents
from su3_levi_equation import active_word, formal, invariants, joint, physical
from su3_levi_kernel import gzip_json
from verify import add_scaled


WEIGHTS = {
    "four_phi1": (3, 4, 3, 1, 1), "three_phi1_phi2": (4, 3, 3, 1, 1),
    "three_phi1_phi3": (4, 4, 2, 1, 1), "three_phi1_minus": (4, 4, 3, 1, 0),
    "three_phi1_plus": (4, 4, 3, 0, 1),
}


def labels(weight, length):
    degree = (*weight[3:], *weight[:3])
    return sorted(bytes(sorted(formal.ID[l] for l in letters)) for letters in contents(degree, length))


def dimension(content):
    character = {0: 1}
    for letter, number in Counter(content).items():
        choose = combinations if formal.ODD[letter] else combinations_with_replacement
        weights = Counter(sum(x) for x in choose((-1, 0, 1), number))
        following = Counter()
        for i, a in character.items():
            for j, b in weights.items():
                following[i+j] += a*b
        character = following
    return character.get(0, 0) - character.get(1, 0)


def frame(weight):
    selected = []
    for content in labels(weight, 7):
        size = dimension(content)
        if not size:
            continue
        rows = [active_word(w) for w in invariants(content)]
        columns = sorted(set().union(*(r.keys() for r in rows)))
        matrix = csr_matrix(np.array([[row.get(k, 0) for k in columns] for row in rows], dtype=np.int64))
        pivots = joint.piv(matrix)
        if len(pivots) != size:
            raise ValueError("two-colour frame misses its Weyl character")
        selected.extend(columns[j] for j in pivots)
    assert len(set(selected)) == len(selected)
    return selected


def polynomial(terms):
    out = {}
    for word, coefficient in terms.items():
        add_scaled(out, active_word(word), coefficient)
    return out


def solve(data, name):
    target = {tuple(bytes.fromhex(t) for t in w): c for w, c in data["coefficients"][name]}
    weight = WEIGHTS[name]
    known_contents = set(labels(weight, 7))
    assert all(bytes(sorted(b"".join(w))) in known_contents for w in target)
    projected = polynomial(target)
    if not projected:
        return dict(channel=name, weight=weight, active_length=7, input_trace_terms=len(target),
                    nonzero_input_monomials=0, denominator=data["denominator"], correction_after_spectator=[],
                    nonzero_residual_monomials=0, kind="polynomial_zero" if target else "formal_zero")
    selected = frame(weight)
    coordinate = {k: i for i, k in enumerate(selected)}
    parents = sorted({w for t in labels(weight, 6) for w in invariants(t)})
    indptr, indices, values = [0], [], []
    for w in parents:
        row = polynomial(formal.q(w))
        entries = sorted((coordinate[k], v) for k, v in row.items() if k in coordinate)
        indices.extend(i for i, v in entries)
        values.extend(v for i, v in entries)
        indptr.append(len(indices))
    matrix = csr_matrix((np.array(values, np.int64), indices, indptr), shape=(len(parents), len(selected)))
    rhs = [-projected.get(k, 0) for k in selected]
    result, solution = physical.certify(matrix, rhs, lambda m: print(name, m, flush=True))
    record = dict(channel=name, weight=weight, active_length=7, parent_trace_products=len(parents),
                  active_coordinates=len(selected), input_trace_terms=len(target),
                  nonzero_input_monomials=len(projected), denominator=data["denominator"], rank=result)
    if result["kind"] == "exact_solution":
        residual = target.copy()
        for i, c in solution.items():
            add_scaled(residual, formal.q(parents[i]), c)
        if polynomial(residual):
            raise ValueError("nonzero polynomial residual at two colours")
        # A central ghost is moved to the left of the active correction.
        sign = -1 if name in ("three_phi1_minus", "three_phi1_plus") else 1
        record["correction_after_spectator"] = [
            [[t.hex() for t in parents[i]], str(sign*c/data["denominator"])] for i,c in sorted(solution.items())]
        record["nonzero_residual_monomials"] = 0
    elif result["kind"] == "exact_separator":
        functional = {selected[i]: c for i, c in solution.items()}
        for w in parents:
            p = polynomial(formal.q(w))
            assert not sum(c*p.get(k, 0) for k, c in functional.items())
        pairing = Fraction(sum(c*projected.get(k, 0) for k, c in functional.items()), data["denominator"])
        assert pairing
        record["functional"] = [[str(k[0]), str(k[1]), str(c)] for k, c in sorted(functional.items())]
        record["pairing_with_Q1"] = str(pairing)
    else:
        raise ValueError(result)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--channel", choices=WEIGHTS, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    data = json.loads(gzip.decompress(args.coefficients.read_bytes()))
    gzip_json(args.out, solve(data, args.channel))
