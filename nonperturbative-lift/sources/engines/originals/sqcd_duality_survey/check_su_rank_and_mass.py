from pathlib import Path
import hashlib
import json
import sympy as sp

nc, nf, r = sp.symbols('Nc Nf r', integer=True)
nt = nf - nc
nlight = nf - r
checks = {}

def zero(name, value):
    result = sp.simplify(value)
    assert result == 0, (name, result)
    checks[name] = str(result)

zero('threshold_exponent', (3*nt-nlight) - ((3*nt-nf)+r))
zero('ADS_branch_number', (nt-nlight) - (r-nc))
zero('low_energy_asymptotic_freedom', (3*nt-nlight) - (2*nt+r-nc))
zero('ADS_mass_dimension', ((3*nt-nf)+r-2*nlight) - 3*(r-nc))
zero('full_rank_condensate_mass_dimension', (3*nt-nf)+nf-3*nt)
zero('electric_mass_vacuum_dimension', (3*nc-nf)+nf-3*nc)
zero('on_shell_W_coefficient', (nt-nf)-(-nc))

phase_cases = []
for nc_value in range(2, 13):
    for nf_value in range(nc_value+2, 3*nc_value):
        nt_value = nf_value-nc_value
        assert (-1)**(nt_value+nf_value) == (-1)**nc_value
        phase_cases.append([nc_value, nf_value])
checks['phase_parity_cases'] = len(phase_cases)

S, L, mu = sp.symbols('S Lambda4 mu', nonzero=True)
mass = sp.symbols('m1:6', nonzero=True)
meson = sp.symbols('M1:6', nonzero=True)
det_m = sp.prod(mass)
det_M = sp.prod(meson)
# Independent determinant substitution for SU(3), Nf=5.
substitution = {Mi:-S/mi for Mi,mi in zip(meson,mass)}
zero('SU3_Nf5_detM', det_M.subs(substitution)+S**5/det_m)
zero('SU3_Nf5_condensate_factorization',
     (S**2-det_M.subs(substitution)/L)
     -S**2*(S**3+L*det_m)/(L*det_m))
zero('SU3_Nf5_Wvac', 2*S+sum(mi*Mi for mi,Mi in zip(mass,meson)).subs(substitution)+3*S)
# Nonzero Lambda4 and all masses make all three roots simple.
poly = sp.Poly(S**3+L*det_m,S)
checks['SU3_Nf5_vacuum_polynomial'] = str(poly.as_expr())
checks['SU3_Nf5_discriminant'] = str(sp.factor(sp.discriminant(poly.as_expr(), S)))
assert sp.factor(sp.discriminant(poly.as_expr(), S)) == -27*(L*det_m)**2

M4 = sp.symbols('det_heavy_M', nonzero=True)
zero('SU3_Nf5_rank4_threshold', (mu**5/L)*M4/mu**4 - mu*M4/L)
checks['SU3_Nf5_rank4_low_energy_exponent'] = 3*2-(5-4)
assert checks['SU3_Nf5_rank4_low_energy_exponent'] == 5

# Exact determinant variation, with generic 2x2 entries, in the derivative order used in the note.
a,b,c,d = sp.symbols('M11 M12 M21 M22')
matrix = sp.Matrix([[a,b],[c,d]])
for row in range(2):
    for col in range(2):
        zero(f'det_variation_{row}{col}', sp.diff(matrix.det(),matrix[row,col])
             - matrix.det()*matrix.inv()[col,row])

root = Path(__file__).resolve().parent
note = root/'su_rank_and_mass.md'
content = note.read_text()
assert sum(line=='$$' for line in content.splitlines()) % 2 == 0
assert not any(token in content for token in ('\\sim','\\approx','\\propto','\\[','\\]'))
checks['note_sha256'] = hashlib.sha256(note.read_bytes()).hexdigest()
(root/'su_rank_and_mass_readback.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
