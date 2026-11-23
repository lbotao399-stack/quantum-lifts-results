"""Check the Taylor constituents and the complete colour frames through degree five."""
from functools import cache
from itertools import combinations, permutations
from pathlib import Path
import gzip
import json

import numpy as np

from b2_schur_complex import multiplicity
from colour import bareiss_residue, matrices
from f4_singlet_characters import decomposition


HERE = Path(__file__).resolve().parent
ZERO = (0, 0, 0, 0)


@cache
def singlet(cycles):
    left, right = [], []
    for length in sorted(cycles, reverse=True):
        (left if sum(left) <= sum(right) else right).append(length)
    a = decomposition("F4", tuple(left))
    b = decomposition("F4", tuple(right))
    return sum(value * b.get(weight, 0) for weight, value in a.items())


def remove_vertical_pair(shape):
    shapes = []
    for i, j in combinations(range(len(shape)), 2):
        short = tuple(v - int(k in (i, j)) for k, v in enumerate(shape))
        if min(short) >= 0 and short == tuple(sorted(short, reverse=True)):
            shapes.append(tuple(v for v in short if v))
    assert len(set(shapes)) == len(shapes)
    return shapes


def low_colour(prime):
    _, adjoint = matrices(HERE)
    form = np.zeros((52, 52), dtype=np.int64)
    data = json.loads(gzip.decompress((HERE / "inputs.json.gz").read_bytes()))
    for i, j, coefficient in data["metric"]:
        form[i, j] = coefficient
    assert np.array_equal(form, form.T)
    assert not np.any((form @ adjoint + adjoint.transpose(0, 2, 1) @ form) % prime)

    def bracket(x, y):
        return (x @ adjoint.reshape(52, -1) % prime).reshape(52, 52) @ y % prime

    def pair(x, y):
        return int((x @ form % prime) @ y % prime)

    def directions(x, degree):
        if degree == 2:
            return [pair(x[0], x[1])]
        if degree == 3:
            return [pair(bracket(x[0], x[1]), x[2])]
        if degree == 4:
            ordinary, marked = [], []
            for i in range(1, 4):
                j, k = sorted(set(range(1, 4)) - {i})
                ordinary.append(pair(x[0], x[i]) * pair(x[j], x[k]) % prime)
                marked.append(pair(bracket(x[0], x[i]), bracket(x[j], x[k])))
            return ordinary + marked[:2]
        result = []
        for i, j in combinations(range(5), 2):
            k, l, m = sorted(set(range(5)) - {i, j})
            result.append(pair(x[i], x[j]) * pair(bracket(x[k], x[l]), x[m]) % prime)
        for i, j, k in permutations((1, 2, 3)):
            result.append(pair(bracket(bracket(bracket(x[0], x[i]), x[j]), x[k]), x[4]))
        return result

    dims = [singlet((1,) * n) for n in range(1, 6)]
    assert dims == [0, 1, 1, 5, 16]
    determinants = []
    for degree in range(2, 6):
        values = []
        for seed in range(dims[degree-1]):
            rng = np.random.Generator(np.random.PCG64(744293 + seed))
            colours = rng.integers(-4, 5, size=(5, 52))
            values.append(directions(colours, degree))
        determinant = bareiss_residue(values, prime)
        assert determinant != 0
        determinants.append(determinant)
    return dict(prime=prime, invariants=dims, ordered_determinants=determinants)


def main():
    target = (4, 4, 3, 3)
    sources = remove_vertical_pair(target)
    assert set(sources) == {(4, 4, 2, 2), (4, 3, 3, 2), (3, 3, 3, 3)}
    rows = []
    for shape in sorted(sources):
        short = shape[2:]
        # In a (2|3) hook, subtracting one from each of the three
        # odd entries removes exactly the top width-three row of the tail.
        target_shape = (shape[0]+1, shape[1]+1, *short[1:]) if short and short[0] == 3 else None
        dims = [multiplicity(shape, length, singlet) for length in (10, 11, 12)]
        n = multiplicity(target_shape, 11, singlet) if target_shape else None
        rows.append(dict(source=list(shape), source_dimensions=dims,
                         first_target=list(target_shape) if target_shape else None,
                         first_target_dimension=n))
    assert rows[0]["source_dimensions"] == [203, 141, 62]
    assert rows[1]["first_target"] == [5, 4, 2] and rows[1]["first_target_dimension"] == 55
    assert [multiplicity((5, 4, 2), length, singlet) for length in (9, 10, 11)] == [834, 311, 55]
    return dict(taylor_products=rows, small_colour=[low_colour(p) for p in (1009, 1013)])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Enable assertions.")
    print(json.dumps(main(), indent=2))
