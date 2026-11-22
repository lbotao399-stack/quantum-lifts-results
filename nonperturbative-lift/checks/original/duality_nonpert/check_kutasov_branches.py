"""Exact specialization of Kutasov--Schwimmer--Seiberg (hep-th/9510222), sec. 6.1.

The script checks the existing duality formula, not a microscopic instanton amplitude.
Internal polynomial variable names are not a new physical normalization.
"""
import json
from pathlib import Path
import sympy as s

root = Path(__file__).resolve().parent
s0 = s.symbols('s0', nonzero=True)
s1, tr_x2, tr_y2 = s.symbols('s1 tr_x2 tr_y2')
nc, nf, magnetic_nc = 5, 4, 3
values = {}
for r1, r2 in [(5, 0), (4, 1), (3, 2)]:
    c1 = s1 / s0 * s.Rational(r2, r1-r2)
    c2 = s1 / s0 * s.Rational(r1, r2-r1)
    assert s.expand(r1*c1+r2*c2) == 0
    assert s.simplify(c1+c2+s1/s0) == 0
    values[f'{r1},{r2}'] = s.factor(r1*c1**2+r2*c2**2)
assert values['5,0'] == 0
assert s.simplify(values['4,1']-s.Rational(20,9)*(s1/s0)**2)==0
assert s.simplify(values['3,2']-30*(s1/s0)**2)==0

# Source equations (6.8), (6.10), and the magnetic characteristic relation.
bar_s0 = -s0
bar_s1 = s.Rational(nc,magnetic_nc)*s1
trace_map = s.Rational(nc,magnetic_nc)*tr_y2 + s.Rational(nc,4)*(s1/s0)**2*(s.Rational(nc**2,magnetic_nc**2)-1)
magnetic_polynomial = tr_y2*(tr_y2-6*(bar_s1/bar_s0)**2)
inverse_map = s.solve(tr_x2-trace_map,tr_y2)[0]
electric_quantum_polynomial = s.expand((tr_x2-values['4,1'])*(tr_x2-values['3,2']))
assert s.simplify(magnetic_polynomial.subs(tr_y2,inverse_map)-s.Rational(9,25)*electric_quantum_polynomial)==0
electric_classical_polynomial = s.expand(tr_x2*electric_quantum_polynomial)
assert s.expand(electric_classical_polynomial.subs(s1,0)-tr_x2**3)==0
assert s.expand(electric_quantum_polynomial.subs(s1,0)-tr_x2**2)==0
assert s.factor(electric_quantum_polynomial.subs(tr_x2,0)) == s.Rational(200,3)*s1**4/s0**4

record={
    'passed':True,
    'source':'hep-th/9510222 sec6.1 eqs6.1,6.5,6.8,6.10,6.12,6.14',
    'Nc':nc,'Nf':nf,'magnetic_Nc':magnetic_nc,
    'trace_values':{k:str(v) for k,v in values.items()},
    'trace_map':str(s.expand(trace_map)),
    'classical_polynomial':str(s.factor(electric_classical_polynomial)),
    'duality_quantum_polynomial':str(s.factor(electric_quantum_polynomial)),
    'removed_branch_value_of_quantum_polynomial':str(s.factor(electric_quantum_polynomial.subs(tr_x2,0))),
    'scope':'Exact specialization of the source duality relation. Adjoint-eigenvalue branches have SQCD moduli; they are not isolated vacua. At s1=0 the nonreduced operator relation uses the source duality map, not a deduction from vacuum values alone. Full electric perturbative survival and a nonperturbative Q primitive are separate gates.'
}
(root/'kutasov_branches_readback.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
