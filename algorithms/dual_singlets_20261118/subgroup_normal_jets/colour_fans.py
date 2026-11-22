"""Resolve the ambient one- and two-loop ordered colour tensors by normal degree."""

from collections import Counter
from pathlib import Path
import json
import sys

from flint import fmpq_mat

ALGORITHMS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ALGORITHMS))
sys.path.insert(0, str(ALGORITHMS / "bmn_primary" / "axis_blocks"))
from bmn_primary_index import G2_ROOTS
from f4_colour import basis as albert_basis
from g2_colours import basis as seven_basis


def commutator(x, y):
    return x * y - y * x


def trace_product(x, y):
    matrix = x * y
    return sum(matrix[i, i] for i in range(matrix.nrows()))


def g2():
    roots, matrices = seven_basis()
    assert len(G2_ROOTS) == 12
    grades = [x % 3 for x, _ in roots]
    assert Counter(grades) == {0: 8, 1: 3, 2: 3}
    # Second simple coroot of the existing seven-dimensional basis.
    return matrices, grades, 3, 2, 1


def f4():
    matrices, _ = albert_basis()
    # Sign of conjugation by diag(1,1,-1) on each Albert coordinate.
    signs = [1] * (3 + 8) + [-1] * 16
    grades = []
    for matrix in matrices:
        rows, columns = matrix.nonzero()
        eigensigns = {signs[i] * signs[j] for i, j in zip(rows, columns)}
        assert len(eigensigns) == 1
        grades.append(int(eigensigns.pop() < 0))
    assert Counter(grades) == {0: 36, 1: 16}
    # Each derivation annihilates the identity, so Tr_27 equals Tr_26.
    return [fmpq_mat(m.tolist()) for m in matrices], grades, 2, 6, grades.index(0)


def colour_fans(name, inputs):
    matrices, grades, period, trace_divisor, probe = inputs
    matrices = list(map(fmpq_mat, matrices))
    size = len(matrices)
    gram = fmpq_mat([[trace_product(a, b) / trace_divisor for b in matrices] for a in matrices])
    inverse = gram.inv()
    duals = []
    for i in range(size):
        dual = fmpq_mat(matrices[0].nrows(), matrices[0].ncols())
        for j, matrix in enumerate(matrices):
            if inverse[i, j]:
                assert (grades[i] + grades[j]) % period == 0
                dual += inverse[i, j] * matrix
        duals.append(dual)
        assert all(trace_product(dual, m) == trace_divisor * int(i == j)
                   for j, m in enumerate(matrices))

    test = matrices[probe]
    assert grades[probe] == 0 and gram[probe, probe]
    left = [commutator(test, t) for t in matrices]
    right = [commutator(t, test) for t in duals]
    first, second = Counter(), Counter()

    def occupancy(*slots):
        return tuple(slots.count(k) for k in range(1, period))

    # Divide by the common external B(test,test) to read the coefficient
    # of the invariant form. The input forms are B(X,Y) and B([X,Y],Z).
    denominator = trace_divisor * gram[probe, probe]
    for i in range(size):
        first[occupancy(grades[i], -grades[i] % period)] += trace_product(left[i], right[i]) / denominator
        for j in range(size):
            middle = commutator(duals[i], matrices[j])
            slots = grades[i], (grades[j] - grades[i]) % period, -grades[j] % period
            second[occupancy(*slots)] += trace_product(commutator(left[i], middle), right[j]) / denominator

    def encode(values):
        return [{"normal_occupations": list(key), "coefficient": str(c)}
                for key, c in sorted(values.items()) if c]

    return {"embedding": name, "trace_divisor": trace_divisor, "graded_dimensions":
            [grades.count(i) for i in range(period)], "external_B": str(gram[probe, probe]),
            "one_loop_colour": encode(first), "two_loop_ordered_colour": encode(second),
            "colour_sums": [str(sum(first.values())), str(sum(second.values()))]}


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run with Python assertions enabled")
    print(json.dumps([colour_fans("A2 in G2", g2()), colour_fans("B4 in F4", f4())], indent=2))
