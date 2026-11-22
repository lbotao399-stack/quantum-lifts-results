"""An integral classical functional on the F4 sixth glueball power."""

import argparse
import gzip
from itertools import combinations, combinations_with_replacement
import json
from math import factorial, gcd
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat

from colour_fans import ALGORITHMS, f4, trace_product
from rank_four_sixth import pfaffian

HERE = Path(__file__).resolve().parent
PAIRS = list(combinations(range(9), 2))


def embedding():
    matrices, grades, _, trace_divisor, _ = f4()
    tangent = [m for m, grade in zip(matrices, grades) if grade == 0]
    # The traceless upper Hermitian 2-by-2 corner of the Albert algebra.
    inclusion = fmpq_mat(27, 9)
    inclusion[0, 0], inclusion[1, 0] = 1, -1
    for i in range(8):
        inclusion[3+i, 1+i] = 1
    coordinates = (inclusion.transpose() * inclusion).inv() * inclusion.transpose()
    vector = [coordinates * t * inclusion for t in tangent]
    assert all(t * inclusion == inclusion * v and v.transpose() == -v
               for t, v in zip(tangent, vector))
    # Invert this matrix to lift any skew 9-by-9 probe into the actual
    # 27-by-27 derivations. The extra invariant line contributes zero trace.
    change = fmpq_mat([[m[i, j] for m in vector] for i, j in PAIRS])
    inverse = change.inv()
    metric = fmpq_mat([[trace_product(a, b) / trace_divisor for b in tangent] for a in tangent])
    for i, j in combinations_with_replacement(range(36), 2):
        assert 2 * metric[i, j] == trace_product(vector[i], vector[j])
        bracket = tangent[i] * tangent[j] - tangent[j] * tangent[i]
        assert bracket * inclusion == inclusion * (vector[i]*vector[j]-vector[j]*vector[i])
    pulled = inverse.transpose() * metric * inverse
    assert pulled == fmpq_mat([[-int(i == j) for j in range(36)] for i in range(36)])
    return np.array(pulled.tolist(), dtype=np.int64)


def target(seed, metric):
    rng = np.random.default_rng(seed)
    coefficients = np.zeros((2, 12, 36), dtype=np.int64)
    for plane in coefficients:
        for colour in plane:
            permutation = rng.permutation(9)
            for i, j in zip(permutation[:8:2], permutation[1:8:2]):
                sign = 2 * int(rng.integers(2)) - 1
                colour[PAIRS.index((min(i, j), max(i, j)))] = sign if i < j else -sign
    mixed = coefficients[0] @ metric @ coefficients[1].T
    return factorial(6) * pfaffian(mixed - mixed.T)


def integer_kernel(matrix):
    hnf, transform = matrix.transpose().hnf(transform=True)
    assert transform * matrix.transpose() == hnf
    rows = [i for i, row in enumerate(hnf.tolist()) if not any(row)]
    return fmpz_mat([transform.tolist()[i] for i in rows])


def compact_dual(kernel, values):
    # Smith divisibility for one column is ordinary integer gcd.
    kernel = kernel.lll()
    hnf, transform = (kernel * values).hnf(transform=True)
    changed = transform * kernel
    assert hnf == changed * values
    particular, *homogeneous = changed.tolist()
    homogeneous = fmpz_mat(homogeneous).lll()
    basis = fmpq_mat(homogeneous)
    projection = (basis * basis.transpose()).solve(basis * fmpq_mat([particular]).transpose())
    nearest = [(2 * int(x.numerator) + int(x.denominator)) // (2 * int(x.denominator))
               for x in projection.entries()]
    result = fmpz_mat([particular]) - fmpz_mat([nearest]) * homogeneous
    assert (result * values)[0,0] == hnf[0,0] > 0
    return list(map(int, result.entries()))


def verify(rebuild=False):
    frame = json.loads(gzip.decompress((ALGORITHMS / "rank_four" / "sixth.json.gz").read_bytes()))
    boundaries = fmpz_mat(frame["boundary_evaluations"])
    metric = embedding()
    native = [target(seed, metric) for seed in frame["seeds"]]
    assert frame["S6_evaluations"] == [64 * value for value in native]
    values = fmpz_mat([native]).transpose()
    kernel = integer_kernel(boundaries)
    assert (boundaries.nrows(), boundaries.ncols(), kernel.nrows()) == (146, 55, 9)
    divisor = gcd(*map(int, (kernel * values).entries()))
    assert divisor == factorial(6)
    path = HERE / "f4_sixth_dual.json"
    if rebuild:
        weights = compact_dual(kernel, values)
        path.write_text(json.dumps({"seeds": frame["seeds"], "weights": weights}, indent=2) + "\n")
    dual = json.loads(path.read_bytes())
    assert dual["seeds"] == frame["seeds"] and all(type(x) is int for x in dual["weights"])
    row = fmpz_mat([dual["weights"]])
    assert (boundaries * row.transpose()).is_zero() and (row * values)[0,0] == divisor
    return dict(group="F4", level=36, spin9_colour_pairs_checked=666,
                rational_boundary_rank=boundaries.ncols()-kernel.nrows(),
                boundary_rows_tested=boundaries.nrows(), probes=len(native),
                integer_image_gcd=divisor, pairing_on_F4_S6=divisor, pairing_on_SO9_S6=64*divisor)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-dual", action="store_true")
    args = parser.parse_args()
    if not __debug__:
        parser.error("Run with assertions enabled")
    print(json.dumps(verify(args.rebuild_dual), indent=2))
