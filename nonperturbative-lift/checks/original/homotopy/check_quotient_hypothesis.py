"""Finite exact matrix check: the operator homotopy does not imply the
ledger's stronger kernel/image contraction without coherent L,h data.

This is a linear-algebra diagnostic for the inference, not a gauge-theory
counterexample. Ordered basis degrees are (-1,0,0,1).
"""
import json
from pathlib import Path
import sympy as s

Qinst = s.zeros(4)
Qinst[2, 0] = 1
Qinst[3, 1] = 1
L = s.zeros(4)
L[2, 1] = 1
h = s.zeros(4)
h[0, 2] = -1
h[1, 3] = -1
assert Qinst**2 == s.zeros(4)
assert L**2 == s.zeros(4)
assert Qinst*L == L*Qinst
assert Qinst*h+h*Qinst == -s.eye(4)
# Exact Q = Lambda^6 Qinst, nu=1 verifies the stated operator assumption.
# ker L / im L is represented by first and fourth basis vectors.
embedding = s.eye(4)[:, [0, 3]]
projection = s.eye(4)[[0, 3], :]
assert projection*Qinst*embedding == s.zeros(2)
assert L*h-h*L != s.zeros(4)
# h sends im L outside im L, preventing descent of the proposed homotopy.
assert (h*L).rank() == 1
assert (L.row_join(h*L)).rank() == 2
report = {
    'Qinst': s.list2numpy(Qinst).tolist(),
    'L': s.list2numpy(L).tolist(),
    'h': s.list2numpy(h).tolist(),
    'kernel_quotient_dimension': 2,
    'induced_Qinst_rank': 0,
    'total_Qinst_cohomology_dimension': 4-2*Qinst.rank(),
    'assumption_scope': 'linear operator identities, no product compatibility asserted',
}
Path('quotient_hypothesis_exact.json').write_text(json.dumps(report, default=int, indent=2))
print(json.dumps(report, default=int))
