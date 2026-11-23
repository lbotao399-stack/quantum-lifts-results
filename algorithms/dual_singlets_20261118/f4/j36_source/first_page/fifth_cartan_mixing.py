"""Weyl dimensions and native lifts of the Cartan target at the fifth glueball charges."""

from functools import cache
from itertools import combinations
import json
from pathlib import Path
from random import Random

from flint import fmpq_mat, fmpz_mat

from source_symmetry import DUAL
from b2_schur_complex import transpose
from f4_cartan_sources import characters, dimensions, weyl_spectra
from single_trace_flag import native26
from trace_representations import cartan


@cache
def schur(shape, traces):
    if len(shape) > 4:
        return 0
    symmetric, _ = characters(traces)
    return int(fmpz_mat([[symmetric[n] if (n := size - i + j) >= 0 else 0
                          for j in range(len(shape))] for i, size in enumerate(shape)]).det())


def count(scalars, spin, marked=False):
    if min((*scalars, *spin)) < 0:
        return 0
    a = tuple(n for n in scalars if n)
    b = transpose(tuple(n for n in spin if n))
    total = sum(size * schur(a, traces) * schur(b, traces) * (traces[0] if marked else 1)
                for traces, size in weyl_spectra().items())
    assert total >= 0 and total % 1152 == 0
    return total // 1152


def test_points():
    weights = cartan(native26())
    def diagonal(vector):
        return [sum(int(w) * x for w, x in zip(weight, vector)) for weight in weights]
    def metric(left, right):
        total = sum(u * v for u, v in zip(left, right))
        assert total % 6 == 0
        return total // 6
    def gram(left, right):
        return fmpz_mat([[metric(a, b) for b in right] for a in left])
    def adjugate(matrix):
        determinant = matrix.det()
        assert determinant
        return (fmpq_mat(matrix).inv() * determinant).numer_denom()[0]

    rng = Random(541036)
    cases = []
    coefficient_rows = [[], []]
    for _ in range(6):
        scalars = [[rng.randrange(-2, 3) for _ in range(4)] for _ in range(3)]
        other = [[rng.randrange(-2, 3) for _ in range(4)] for _ in range(3)]
        x, z = list(map(diagonal, scalars)), list(map(diagonal, other))
        xx = gram(x, x)
        if not xx.det() or not gram(z, z).det():
            continue
        adj = adjugate(xx)
        assert xx * adj == fmpz_mat([[xx.det() * (i == j) for j in range(3)] for i in range(3)])
        mixed = gram(x, z[:1])
        sixth = 0
        for u, v, w, t in zip(*x, z[0]):
            col = fmpz_mat([[u], [v], [w]])
            sixth += t**2 * (col.transpose() * adj * col)[0, 0] ** 2
        one = [int(xx.det()**2 * metric(z[0], z[0])),
               int(xx.det() * (mixed.transpose() * adj * mixed)[0, 0]), int(sixth)]
        three = [int(xx.det() * gram(z, z).det()), int(gram(x, z).det()**2)]
        cases.append(dict(scalars=scalars, dotted_test_vectors=other,
                          one_dotted_pair=one, three_dotted_pairs=three))
        for rows, values in zip(coefficient_rows, (one, three)):
            rows.append(values)
    witnesses = []
    for rows, rank in zip(coefficient_rows, (3, 2)):
        for selected in combinations(range(len(rows)), rank):
            determinant = int(fmpz_mat([rows[i] for i in selected]).det())
            if determinant:
                witnesses.append(dict(rows=selected, determinant=str(determinant)))
                break
        else:
            raise AssertionError("native Cartan traces failed to span")
    return dict(evaluations=cases, witnesses=witnesses)


def calculate():
    archived = json.loads((DUAL / "f4/cartan.json").read_bytes())
    old = {row["t"]: sum(row["dimensions"].values()) for row in archived["cartan_dimensions"] if row["power"] == 5}
    rows = []
    for k in range(6):
        p = 5 - k
        target = count((p, p, p), (k, k))
        source = [count((p, p - 1, p - 1), (k, k), True),
                  count((p, p, p - 1), (k, k - 1), True),
                  count((p, p, p), (k, k - 2), True)]
        weight = (k, k, p, p, p)
        assert sum(source) == old[k]
        assert target == sum(dimensions(weight, sum(weight)).values())
        rows.append(dict(dotted_pairs=k, scalars_per_flavour=p, target=target,
                         source_by_second_jet=source, source_total=sum(source)))
    assert [r["target"] for r in rows] == [0, 3, 0, 2, 0, 0]
    return dict(weyl_order=sum(weyl_spectra().values()), target_twice_energy=30, target_twice_R=[10]*3,
                source_twice_energy=29, source_twice_R=[9]*3, cartan_dimensions=rows, native=test_points())


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(), indent=2))
