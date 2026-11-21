"""Bound the kernel of the six returned odd Ward actions on the triality doublet."""
import argparse
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

import mixed_incoming as m
from so8_level24_descendant import check_letter, substitute
from so8_triality import representations
from so8_triality_page import all_rows, probe, rank
from trace_descent import dimensions


@lru_cache(maxsize=40000)
def returned_word(word):
    total = Counter()
    for spin in range(2):
        for flavour in range(3):
            first = substitute({word: 1}, spin, flavour, False)
            if first:
                total.update(substitute(first, spin, flavour, True))
    return tuple((w, n) for w, n in total.items() if n)


def returned_row(row):
    result = Counter()
    for word, coefficient in row.items():
        for w, factor in returned_word(word):
            result[w] += coefficient * factor
    return {w: n for w, n in result.items() if n}


def run(grade, deficit, seed):
    weight = (grade, grade) + (4 - grade,) * 3
    length = sum(weight) - deficit
    parents = all_rows(weight, length - 1)
    centre = all_rows(weight, length)
    alphabet = {letter for row in centre for word in row for trace in word for letter in trace}
    for letter in alphabet:
        for spin in range(2):
            for flavour in range(3):
                check_letter(letter, spin, flavour, False)
                check_letter(letter, spin, flavour, True)
    boundaries = [m.q(row) for row in parents]
    outgoing = [m.q(row) for row in centre]
    returns = [returned_row(row) if grade else {} for row in centre]
    assert all(not m.q(row) for row in boundaries)
    print("grade", grade, "length", length, "formal rows", len(parents), len(centre),
          "return terms", sum(map(len, returns)), flush=True)
    dims = dimensions(weight, length, 8, "SO")
    next_dims = dimensions(weight, length + 1, 8, "SO")
    spinors = representations()
    together = probe(centre + boundaries + returns, dims, seed, spinors)
    following = probe(outgoing, next_dims, seed + 17417, spinors)[2]
    n, b = len(centre), len(boundaries)
    full, invariant, standard = together
    cochains = rank(standard[:n])
    invariant_dimension = rank(invariant[:n])
    assert rank(full[:n]) == cochains + invariant_dimension == sum(dims.values())
    incoming = standard[n:n+b]
    action = standard[n+b:]
    matrix = np.block([[following, action],
                       [np.zeros((b, following.shape[1]), dtype=np.int64), incoming]])
    combined = rank(matrix)
    before, after = rank(incoming), rank(following)
    assert before + after <= combined <= cochains
    return dict(grade=grade, weight=weight, length=length, deficit=deficit, prime=m.P, seed=seed,
                letter_types=len(alphabet), source_rows=n, predecessor_rows=b,
                target_character=dims, following_character=next_dims,
                invariant_dimension=invariant_dimension, standard_even_dimension=cochains,
                incoming_rank_lower=before, outgoing_rank_lower=after, joint_rank_lower=combined,
                classical_upper=cochains-before-after, ward_kernel_upper=cochains-combined)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grade", type=int, choices=range(5), required=True)
    p.add_argument("--deficit", type=int, choices=(1, 3), required=True)
    p.add_argument("--seed", type=int, default=873731)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    data = run(args.grade, args.deficit, args.seed)
    m.write(args.out, data)
    print(data, flush=True)
