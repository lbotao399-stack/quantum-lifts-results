"""Highest weights in the terminal N=4 super-commuting ring for USp(2n)."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import json
from math import factorial, prod
from pathlib import Path


@lru_cache(None)
def determinants(rank):
    """Multiplicity of each det(1-x w) for signed permutations on the Cartan."""
    counts = Counter()
    for perm in permutations(range(rank)):
        for signs in product((1, -1), repeat=rank):
            polynomial = [1]
            seen = set()
            for start in range(rank):
                if start in seen:
                    continue
                position, length, sign = start, 0, 1
                while position not in seen:
                    seen.add(position)
                    length += 1
                    sign *= signs[position]
                    position = perm[position]
                polynomial += [0] * length
                for power in reversed(range(length, len(polynomial))):
                    polynomial[power] -= sign * polynomial[power-length]
            counts[tuple(polynomial)] += 1
    assert sum(counts.values()) == 2**rank * factorial(rank)
    return counts


@lru_cache(None)
def symmetric(determinant, degree):
    if degree < 0:
        return 0
    if degree == 0:
        return 1
    return -sum(c * symmetric(determinant, degree-i)
                for i, c in enumerate(determinant[1:], 1) if i <= degree)


def exterior(determinant, degree):
    return (-1)**degree * determinant[degree] if 0 <= degree < len(determinant) else 0


@lru_cache(None)
def weight_dimension(rank, spin, flavours):
    classes = determinants(rank)
    total = sum(count * prod(exterior(d, n) for n in spin) *
                prod(symmetric(d, n) for n in flavours) for d, count in classes.items())
    dim, rem = divmod(total, sum(classes.values()))
    assert rem == 0 and dim >= 0
    return dim


def highest(rank, spin, flavours):
    result = 0
    for perm in permutations(range(3)):
        sign = (-1)**sum(i > j for k, i in enumerate(perm) for j in perm[k+1:])
        shifted = tuple(n+i-k for k, (n, i) in enumerate(zip(flavours, perm)))
        result += sign * (weight_dimension(rank, spin, shifted) -
                          weight_dimension(rank, (spin[0]+1, spin[1]-1), shifted))
    assert result >= 0
    return result


def column(rank, weight):
    return [dict(derivative_degree=t, length=sum(weight)-t,
                 highest_weight_dimension=highest(rank, (t, t), tuple(n-t for n in weight)))
            for t in range(min(weight)+1)]


def run(public_table):
    public = json.loads(public_table.read_bytes())
    assert public["jmax"] == 28 and public["checks_failed"] == 0
    blocks = {tuple(row["lam"]): row for row in public["columns"] if row["twoJ"] == 0}
    tested = 0
    for a in range(1, 15):
        for b in range(a+1):
            for c in range(b+1):
                if a+b+c > 14:
                    continue
                weight = a, b, c
                h = highest(2, (0, 0), weight)
                assert h == blocks.get(weight, {}).get("h", 0)
                tested += 1
    return dict(
        public_terminal_weights_checked=tested,
        rank_two=[dict(weight=w, terminal_grades=column(2, w)) for w in
                  ((10,4,4), (11,4,4), (12,4,4), (10,5,5), (13,4,4), (11,5,5))],
        rank_three_omega=column(3, (8,3,3)),
        weyl_group_orders={r: sum(determinants(r).values()) for r in (2,3)},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--so5-table", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(run(args.so5_table), indent=2) + "\n")
