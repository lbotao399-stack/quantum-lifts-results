#!/usr/bin/env python3
"""Exact complete USp(4) Konishi source/target and physical mixing audit.

This is a channel calculation, not an all-order physical comparison theorem.
The quantum normalization is fixed by the physical multiplet calculation in
README.md; the independent short-distance integral checks its magnitude.
"""
from pathlib import Path
import sys
CHANNEL_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CHANNEL_DIR.parent))
from independent_checks import *
from fortuity.full.q1 import pair_rule
import sympy as sp

def complete_words(weight, length):
    words = set()
    for content in contents(weight, length):
        for word in multitraces_with_content(content):
            sign, value = canonical(word)
            if sign:
                words.add(value)
    return sorted(words)

def covariant_source(flavour):
    complement = tuple(i for i in range(3) if i != flavour)
    indices = complement + (flavour,)
    orientation = (-1)**sum(x > y for i, x in enumerate(indices)
                             for y in indices[i + 1:])
    sign, word = canonical((((0, 0, (flavour,)), (0, 0, complement)),))
    return {word: sign * orientation}

def native_quantum_bilinear(row, adjoint_casimir):
    result = Counter()
    for word, value in row.items():
        assert len(word) == 1 and len(word[0]) == 2
        left, right = word[0]
        for coefficient, a, b in pair_rule(left[2], right[2]):
            # Exact compact-matrix identity, independently checked below:
            # sum_C tr([T_C,a][T_C,b])=-C_adj tr(a b).
            sign, target = canonical(((a, b),))
            if sign:
                result[target] -= value * coefficient * adjoint_casimir * sign
    return {word: value for word, value in result.items() if value}

def compact_colour_check():
    raw, _ = generators()
    matrices = [sp.Matrix(a.tolist()) * sp.I**phase / sp.sqrt(norm)
                for a, phase, norm in raw]
    pairings = []
    for left in matrices:
        row = []
        for right in matrices:
            value = sum(sp.trace((a * left - left * a) *
                                (a * right - right * a)) for a in matrices)
            row.append(sp.simplify(value))
        pairings.append(row)
    assert sp.Matrix(pairings) == -6 * sp.eye(10)
    return [[str(value) for value in row] for row in pairings]

def physical_sector_weights():
    # Exact solutions of all five physical charge constraints, including
    # candidate one-letter traces which vanish by the simple gauge algebra.
    result = {}
    for name, twice_energy, twice_R, twice_left in (
            ('source', 5, (1, 1, 1), 1),
            ('target', 6, (2, 2, 2), 0)):
        solutions = []
        for length in range(1, twice_energy + 1):
            for a, b in product(range(twice_energy + 1), repeat=2):
                for theta in product(range(twice_energy + 1), repeat=3):
                    total = sum(theta)
                    if (2 * (a + b) + total + length != twice_energy
                            or a != b or a + b + total - length != twice_left
                            or tuple(length - total + 2 * value for value in theta) != twice_R):
                        continue
                    weight = (a, b) + theta
                    words = complete_words(weight, length)
                    solutions.append({'native_weight': list(weight),
                                      'length': length,
                                      'independent_trace_words': len(words)})
        result[name] = solutions
    assert result == {
        'source': [{'native_weight': [1, 1, 0, 0, 0], 'length': 1,
                    'independent_trace_words': 0},
                   {'native_weight': [0, 0, 1, 1, 1], 'length': 2,
                    'independent_trace_words': 3}],
        'target': [{'native_weight': [1, 1, 0, 0, 0], 'length': 2,
                    'independent_trace_words': 1},
                   {'native_weight': [0, 0, 1, 1, 1], 'length': 3,
                    'independent_trace_words': 1}]}
    return result

def spinor_check():
    # Preserve covariant spinor contractions in the proof; this independent
    # exact check verifies every spinor entry for all homogeneous monomials
    # of the separation vector, without an angular average.
    sigma = [sp.eye(2), sp.Matrix([[0, sp.I], [sp.I, 0]]),
             sp.Matrix([[0, 1], [-1, 0]]), sp.diag(sp.I, -sp.I)]
    epsilon = sp.Matrix([[0, 1], [-1, 0]])
    coefficients = []
    for mu in range(4):
        for nu in range(mu, 4):
            value = sigma[mu] * epsilon * sigma[nu].T
            if mu != nu:
                value += sigma[nu] * epsilon * sigma[mu].T
            expected = epsilon if mu == nu else sp.zeros(2)
            assert sp.simplify(value - expected) == sp.zeros(2)
            coefficients.append({'mu': mu, 'nu': nu,
                                 'spinor_tensor': [[str(x) for x in value.row(i)]
                                                   for i in range(2)]})
    return coefficients

def main():
    started = time.monotonic()
    sector_weights = physical_sector_weights()
    source = [covariant_source(i) for i in range(3)]
    scalar = {complete_words((0, 0, 1, 1, 1), 3)[0]: 1}
    gaugino = {complete_words((1, 1, 0, 0, 0), 2)[0]: 1}
    target = [scalar, gaugino]
    assert len({next(iter(row)) for row in source}) == 3
    for row in source:
        assert q(row) == {word: 2 * value for word, value in scalar.items()}
    assert not q(scalar) and not q(gaugino)
    colour = compact_colour_check()
    for row in source:
        assert native_quantum_bilinear(row, 6) == {
            word: -12 * value for word, value in gaugino.items()}
    source_gram, _ = gram(source, label='F03 source')
    target_gram, _ = gram(target, label='F03 complete target')
    assert source_gram == 10 * fmpq_mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert target_gram == fmpq_mat([[15, 0], [0, 10]])
    # Entire independent colour-coordinate Wick calculation on all five rows.
    raw, _ = generators()
    direct_rows = []
    for row in source + target:
        result = Counter()
        for word, value in row.items():
            for term, coefficient in coefficient_polynomial(word, raw).items():
                result[term] += value * coefficient
        direct_rows.append(dict(result))
    direct = [[str(direct_pair(left, right, raw)) for right in direct_rows]
              for left in direct_rows]
    assert direct == [['10', '0', '0', '0', '0'], ['0', '10', '0', '0', '0'],
                      ['0', '0', '10', '0', '0'], ['0', '0', '0', '15', '0'],
                      ['0', '0', '0', '0', '10']]
    classical = sp.Matrix([[2, 2, 2], [0, 0, 0]])
    quantum = sp.Matrix([[0, 0, 0], [-12, -12, -12]])
    source_metric = 10 * sp.eye(3)
    target_metric = sp.diag(15, 10)
    radial_classical = source_metric.inv() * classical.T * target_metric
    assert radial_classical == sp.Matrix([[3, 0], [3, 0], [3, 0]])
    assert (radial_classical * classical).eigenvals() == {0: 2, 18: 1}
    kernel = sp.Matrix([[1, 0], [-1, 1], [0, -1]])
    assert classical * kernel == sp.zeros(2)
    assert quantum * kernel == sp.zeros(2)
    g = sp.Symbol('g_YM', real=True)
    physical = g * classical / (4 * sp.pi) - sp.I * g**2 * quantum / (32 * sp.pi**2)
    physical_adjoint = source_metric.inv() * physical.conjugate().T * target_metric
    protected = sp.Matrix([sp.I * g / (2 * sp.pi), 1])
    long = sp.Matrix([1, 3 * sp.I * g / (4 * sp.pi)])
    orthogonality = (protected.conjugate().T * target_metric * long)[0]
    assert sp.simplify(orthogonality) == 0
    assert physical == sp.Matrix.hstack(*(g * long / (2 * sp.pi) for _ in range(3)))
    assert sp.simplify(physical_adjoint * protected) == sp.zeros(3, 1)
    assert physical * kernel == sp.zeros(2)
    # Exact Feynman-parameter evaluation of the finite Wilson-line/Yukawa
    # residue. The radial integral is pi^2/(2*Delta), evaluated before the
    # parameter cancellation; no numerical limiting extrapolation is used.
    t = sp.Symbol('t', positive=True)
    separation_squared = sp.Symbol('separation_squared', positive=True)
    radial = sp.pi**2 / (2 * t * (1 - t) * separation_squared)
    convolution = sp.integrate(sp.cancel(2 * t * (1 - t) * radial), (t, 0, 1))
    assert convolution == sp.pi**2 / separation_squared
    # Magnitudes from CFLT action (D.9), Wilson line, and two propagators.
    loop_magnitude = (2 * 128 / g**2 * (g**2 / (64 * sp.pi**2))**2
                      * sp.pi**2 * sp.Rational(1, 2) * 6)
    assert sp.simplify(loop_magnitude) == 3 * g**2 / (16 * sp.pi**2)
    spinors = spinor_check()
    record = {
        'passed': True,
        'scope': 'Complete channel algebra and physical first-quantum-order mixing checks; acceptance of local-product/current comparison requires the accompanying proof.',
        'gauge_group': 'USp(4)',
        'source_charges': {'E': '5/2', 'R': ['1/2', '1/2', '1/2'], 'J_L': '1/2', 'J_R': '0'},
        'target_charges': {'E': '3', 'R': ['1', '1', '1'], 'J_L': '0', 'J_R': '0'},
        'complete_native_weight_exhaustion': sector_weights,
        'source_native_trace_polynomials': export_rows(source),
        'target_native_trace_polynomials': export_rows(target),
        'source_gram': export_matrix(source_gram),
        'target_gram': export_matrix(target_gram),
        'all_finite_rank_null_relations_after_cyclic_and_symplectic_reversal': [],
        'native_classical_matrix': [[str(x) for x in classical.row(i)] for i in range(2)],
        'native_pair_rule_matrix_before_physical_normalization': [[str(x) for x in quantum.row(i)] for i in range(2)],
        'exact_source_kernel_basis': [[str(x) for x in kernel.row(i)] for i in range(3)],
        'physical_supercharge_coefficients_through_g_squared': [[str(x) for x in physical.row(i)] for i in range(2)],
        'physical_radial_adjoint_coefficients_through_g_squared': [[str(x) for x in physical_adjoint.row(i)] for i in range(3)],
        'protected_target_direction_through_g': [str(x) for x in protected],
        'long_target_direction_through_g': [str(x) for x in long],
        'protected_long_pairing_at_the_required_order': str(orthogonality),
        'Konishi_one_loop_anomalous_dimension': '9*g_YM**2/(4*pi**2)',
        'independent_compact_colour_matrix': colour,
        'independent_all_25_Wick_pairings': direct,
        'independent_spinor_convolution_coefficients': spinors,
        'exact_point_split_convolution_vector_coefficient': str(convolution),
        'point_split_anomaly_magnitude_in_CFLT_field_convention': str(sp.simplify(loop_magnitude)),
        'seconds': time.monotonic() - started,
    }
    write(CHANNEL_DIR / 'channel_receipt.json', record)
    print(json.dumps({'passed': True, 'source_dimension': 3, 'target_dimension': 2,
                      'source_null_relations': 0, 'target_null_relations': 0,
                      'independent_Wick_pairings': 25, 'seconds': record['seconds']}))

if __name__ == '__main__':
    main()
