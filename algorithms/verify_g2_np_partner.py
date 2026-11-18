"""Replay the rational G2 partner ranks and the six-product polynomial identity."""
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from sympy import Matrix, zeros


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "su3_independent_prime_20261118" / "runtime" / "work" / "sp3-native.nosync"))
from sp_environment import canonical, fw


def matrix(rows):
    return Matrix([[Fraction(value) for value in row] for row in rows])


def polynomial(rows):
    result = {}
    for coefficient, word in rows:
        traces = tuple(tuple((a, b, tuple(theta)) for a, b, theta in trace) for trace in word)
        result[traces] = Fraction(coefficient)
    return result


def run():
    path = HERE.parent / "data" / "N4_G2_np_partner_replay_20261118.json.gz"
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "a3d95d83b6e8ea512a972730f1747dcfa172a47068ff0dbc4b11e086bb7d7e5a"
    data = json.loads(gzip.decompress(raw))

    for key, dimension, boundaries in (("S2_evaluations", 2, 1), ("Gram_evaluations", 4, 3)):
        assert len(data[key]) == 2
        for saved_rows in data[key]:
            rows = matrix(saved_rows)
            assert rows.rows == dimension + boundaries + 1
            assert rows[:dimension, :].rank() == dimension
            assert rows[dimension:-1, :].rank() == boundaries
            assert rows[dimension:, :].rank() == boundaries + 1

    source = {k: matrix(v) for k, v in data["source"].items()}
    coefficients = matrix([data["source_coefficients"]])
    frame = source["source_values"]
    boundaries = source["boundaries"]
    assert frame.shape == (5, 5) and frame.rank() == 5
    assert boundaries.rank() == 1 and source["Q0"].rank() == 3
    assert boundaries * frame.inv() * source["Q0"] == zeros(boundaries.rows, 3)
    assert coefficients * source["Q0"] == zeros(1, 3)
    assert boundaries.col_join(coefficients * frame).rank() == 2
    coefficient = Fraction(data["correction_coefficient"])
    assert source["Q1"] + coefficient * source["Q0_correction"] == zeros(1, 30)

    combined = defaultdict(Fraction)
    assert len(data["six_products"]) == 6
    for pair in data["six_products"]:
        left, right = polynomial(pair["first"]), polynomial(pair["second"])
        assert not canonical(fw.Q_expr(left)) and not canonical(fw.Q_expr(right))
        for first_word, first_coefficient in left.items():
            for second_word, second_coefficient in right.items():
                combined[first_word + second_word] += first_coefficient * second_coefficient
    source_polynomial = polynomial(data["representative"])
    assert len(source_polynomial) == 18 and canonical(combined) == source_polynomial
    return {"source_cohomology_dimension": 5 - 1 - 3, "independent_target_classes_checked": 2,
            "first_correction_coefficient": str(coefficient), "rational_residual_coordinates": 30,
            "closed_bilinear_factor_occurrences": 12, "reconstructed_trace_terms": 18,
            "nonperturbative_matrix_computed": False}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
