"""Exact global Jacobian-ideal identity in Seiberg's constrained effective model."""
import json
from pathlib import Path
import sympy as s

X, M11, M12, M21, M22, B, Bt, scale, mass = s.symbols(
    "X M11 M12 M21 M22 B Bt Lambda4 m2"
)
W = X*(M11*M22-M12*M21-B*Bt-scale) + mass*M22
combination = (
    scale*s.diff(W, M22)
    - M11**2*s.diff(W, M11)
    - M11*M12*s.diff(W, M12)
    - M11*B*s.diff(W, B)
    + M11*X*s.diff(W, X)
)
residual = s.expand(combination-mass*scale)
assert residual == 0
# At zero instanton factor, the critical locus includes M11=1, X=-m2.
zero_scale_point = {scale: 0, M11: 1, X: -mass,
                    M12: 0, M21: 0, M22: 0, B: 0, Bt: 0}
critical_residuals = [s.expand(s.diff(W,v).subs(zero_scale_point))
                      for v in [X,M11,M12,M21,M22,B,Bt]]
assert critical_residuals == [0]*7
out = {"jacobian_identity_residual": str(residual),
       "constant_in_jacobian_ideal": "m2 * Lambda^4",
       "zero_scale_critical_point_residuals": list(map(str,critical_residuals)),
       "scope": "global constrained effective Koszul model; no inverse meson",
       "status": "PASS"}
Path(__file__).with_name("global_ads_identity_readback.json").write_text(
    json.dumps(out, indent=2)+"\n"
)
print(json.dumps(out))
