"""Native root centralizers and scalar-determinant controls for the F4 trace test."""
from collections import Counter
import json
from math import factorial
from random import Random

from flint import fmpq_mat, fmpz_mat
import numpy as np

from colour import Colour
from scalar_plane_ward import blocks, coordinate_groups, plane
from single_trace_flag import native26


def determinant(a, b, c):
    return (a[0] * (b[1]*c[2] - b[2]*c[1]) + a[1] * (b[2]*c[0] - b[0]*c[2])
            + a[2] * (b[0]*c[1] - b[1]*c[0]))


def check_ward(d, power, modulus):
    seed = 361144
    coefficients, degrees, grades, _ = plane(d, modulus, seed)
    rows, cols = coordinate_groups(d, power, modulus, degrees, grades)
    rng = Random(6011)
    scalars = [[rng.randint(-3, 3) for _ in range(3)] for _ in range(d)]
    ghost = [rng.randint(-3, 3) for _ in range(52)]
    t = native26().reshape(52, -1)
    abc = (np.array(scalars).T @ np.array(coefficients) @ t).reshape(3, 26, 26)
    u = (np.array(ghost) @ t).reshape(26, 26)
    gram = fmpq_mat((t @ t.reshape(52, 26, 26).transpose(0, 2, 1).reshape(52, -1).T).tolist())
    comm = (abc @ u - u @ abc).transpose(0, 2, 1).reshape(3, -1)
    recovered = gram.solve(fmpq_mat((t @ comm.T).tolist()))
    assert all(x.denominator == 1 for x in recovered.entries())
    raw = [[int(x) for x in r] for r in recovered.tolist()]
    def value(col):
        return determinant(*(scalars[i] for i in col)) if col else 1
    count = 0
    for w, h, width, entries in blocks(d, power, modulus, seed):
        v = [value(s[:3]) * value(s[3:]) * ghost[i] for s, i in cols[w]]
        observed = [0]*h
        for i, j, x in entries:
            observed[i] += x*v[j]
        expected = [value(s) * determinant(scalars[i], scalars[j], raw[k]) for s, (i,j), k in rows[w]]
        assert len(v) == width and observed == expected
        count += h
    return dict(colour_plane=d, determinant_degree=power, native_rows=count)


def calculate():
    c = Colour()
    lengths = [int(c.scalar_product(w, w)) for w in c.weights]
    records = Counter()
    examples = {}
    for e in range(4, 52):
        f = c.opposite[e]
        roots = [v for v in range(4, 52) if not c.br[e][v] and not c.br[f][v]]
        short = sum(lengths[i] == 1 for i in roots)
        long = sum(lengths[i] == 2 for i in roots)
        records[lengths[e], short, long] += 1
        examples.setdefault(lengths[e], (e, roots))
    assert records == {(2, 12, 6): 24, (1, 0, 12): 24}
    torus = native26()[:4]
    gram = fmpz_mat([[int(np.trace(x @ y)) for y in torus] for x in torus])
    # Order all four U before all four V as in the native exterior coordinates.
    determinant = int(gram.det())
    scalar_minor = int(fmpz_mat([r[:3] for r in gram.tolist()[:3]]).det())
    value = (-1)**6 * factorial(4) * determinant * scalar_minor
    assert value
    return dict(root_sl2_centralizers=[
        dict(root_norm_squared=k, count=v, short_root_spaces=short, long_root_spaces=long,
             sample_root=examples[k][0], sample_centralizing_roots=examples[k][1])
        for (k,short,long),v in sorted(records.items())],
        cartan_defining_trace_gram=[[int(x) for x in r] for r in gram.tolist()],
        cartan_scalar_gram_determinant=scalar_minor, cartan_glueball_gram_determinant=determinant,
        cartan_glueball_fourth_power_times_scalar_gram=value,
        full_defining_matrix_ward_checks=[check_ward(4,1,1),check_ward(5,1,1),
                                          check_ward(6,2,7),check_ward(7,2,7)])


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(), indent=2))
