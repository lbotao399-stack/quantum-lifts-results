"""Exact, bounded check of the one-matrix SU(5) invariant quotient."""
import json
from pathlib import Path

import sympy as sp

out = Path(__file__).with_name("su5_invariant_quotient_readback.json")
s2, s3, s4, s5 = sp.symbols("s2 s3 s4 s5")
# s_j is only a machine identifier for Tr(X**j).
power_sums = {1: sp.Integer(0), 2: s2, 3: s3, 4: s4, 5: s5}
elementary = {0: sp.Integer(1)}
for degree in range(1, 6):
    elementary[degree] = sp.expand(sum(
        (-1) ** (k - 1) * elementary[degree-k] * power_sums[k]
        for k in range(1, degree + 1)
    ) / degree)

# Newton recurrence for p6 in dimension five, with p1=0.
s6 = sp.expand(sum(
    (-1) ** (k + 1) * elementary[k] * power_sums[6-k]
    for k in range(1, 6)
))
generators = [s3, s4-s2**2/5, s5-s2*s3/5, s6-s2*s4/5]
basis = sp.groebner(generators, s5, s4, s3, s2, order="lex", domain=sp.QQ)
expected = [s5, s4-s2**2/5, s3, s2**3]
assert list(basis) == expected
assert basis.reduce(s2**2)[1] == s2**2
assert basis.reduce(s2**3)[1] == 0

# Explicit generator combination producing the cubic, independent of Groebner.
cubic_combination = sp.expand(
    generators[3] - sp.Rational(11, 20)*s2*generators[1]
    - sp.Rational(1, 3)*s3*generators[0]
)
assert cubic_combination == -sp.Rational(3, 200)*s2**3

result = {
    "coefficient_domain": "QQ",
    "variable_dictionary": {str(p): f"Tr(X^{n})" for n,p in power_sums.items() if n > 1},
    "elementary_symmetric_polynomials": {str(k): str(v) for k,v in elementary.items()},
    "Tr_X6": str(s6),
    "trace_F_generators_j_1_to_4": list(map(str, generators)),
    "lex_order": ["s5", "s4", "s3", "s2"],
    "reduced_groebner_basis": list(map(str, basis)),
    "square_normal_form": str(basis.reduce(s2**2)[1]),
    "cube_normal_form": str(basis.reduce(s2**3)[1]),
    "cubic_generator_combination": str(cubic_combination),
    "conclusion": "C[Tr(X^2)] / ((Tr(X^2))^3); square nonzero",
    "scope": "Completeness of trace generators requires the Reynolds/covariant proof in su5_invariant_quotient.md.",
}
out.write_text(json.dumps(result, indent=2) + "\n")
assert json.loads(out.read_text()) == result
print(json.dumps(result, indent=2))
