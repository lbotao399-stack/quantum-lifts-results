"""Relative SU(2) bilinear jets and the half-BPS scalar descendant kernel."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial

import sympy as sp


def sign(items):
    return (-1)**sum(left > right for left, right in combinations(items, 2))


def add_exterior(polynomial, fields, coefficient):
    if len(set(fields)) == len(fields):
        polynomial[tuple(sorted(fields))] += coefficient * sign(fields)


@lru_cache(maxsize=None)
def differential(field):
    species, order, color = field
    result = Counter()
    for first, second in permutations([i for i in range(3) if i != color]):
        epsilon = sign((color, first, second))
        for p in range(1, order+1):
            q = order - p
            if species == "b":
                add_exterior(result, (("c", p, first), ("b", q, second)), epsilon)
            elif q > 0:
                add_exterior(result, (("c", p, first), ("c", q, second)), epsilon * sp.S.Half)
    return result


def boundary(polynomial):
    result = Counter()
    for fields, coefficient in polynomial.items():
        for j, field in enumerate(fields):
            for wedge, value in differential(field).items():
                add_exterior(result, fields[:j] + wedge + fields[j+1:],
                             coefficient * (-1)**j * value)
    return {fields: v for fields, v in result.items() if v}


def bilinears(spin):
    # Modes are the Taylor coefficients in the source. A b-mode has one
    # additional unit of spin, a c-mode starts at derivative order one.
    modes = [(species, order) for species in ("b", "c")
             for order in range(int(species == "c"), spin+1)]
    for left, right in combinations(modes, 2):
        if left[1] + right[1] + (left[0] == "b") + (right[0] == "b") != spin:
            continue
        polynomial = Counter()
        for color in range(3):
            add_exterior(polynomial, ((*left, color), (*right, color)), 1)
        yield left, right, boundary(polynomial)


def pure_cubic(order):
    result = Counter()
    for p in range(order+1):
        for q in range(order-p+1):
            r = order - p - q
            for i, j, k in permutations(range(3)):
                add_exterior(result, (("b", p, i), ("b", q, j), ("b", r, k)), sign((i,j,k)))
    return {fields: v for fields, v in result.items() if v}


def scalar_descendants():
    length, first, second, along_first, along_second = sp.symbols("L U V a d")
    for degree in (2, 3, 4):
        # (x-au-dv)^2 with a totally isotropic pair of translation vectors.
        correlator = (length - 2*along_first*first - 2*along_second*second)**(-degree)
        for order in range(7):
            for r in range(order+1):
                value = correlator.diff(along_first, r, along_second, order-r).subs(
                    {along_first: 0, along_second: 0})
                coefficient = 2**order * factorial(degree+order-1) // factorial(degree-1)
                expected = coefficient*first**r*second**(order-r)/length**(degree+order)
                assert sp.cancel(value - expected) == 0
    print("N=4 two-point derivative coefficients verified at scalar degrees 2,3,4 and orders 0 to 6.")


def main():
    for spin in range(2, 17):
        columns = list(bilinears(spin))
        rows = sorted(set().union(*(image.keys() for _, _, image in columns)))
        matrix = sp.SparseMatrix(len(rows), len(columns),
                                 {(j, k): image[row] for k, (_, _, image) in enumerate(columns)
                                  for j, row in enumerate(rows) if row in image})
        jet = sp.Matrix([right[1] if (left[0], right[0]) == ("b", "c") else 0
                         for left, right, _ in columns])
        assert matrix * jet == sp.zeros(len(rows), 1) and any(jet)
        assert matrix.rank() == len(columns) - 1
        assert all(any(field[0] == "c" for field in row) for row in rows)
        if spin >= 3:
            cubic = pure_cubic(spin-3)
            assert cubic and not boundary(cubic)
    print("Full relative SU(2) bilinear kernels have dimension one at spins 2 through 16.")
    print("Each is the derivative of X. The pure-b derivatives of Y are closed and nonzero.")
    scalar_descendants()


if __name__ == "__main__":
    main()
