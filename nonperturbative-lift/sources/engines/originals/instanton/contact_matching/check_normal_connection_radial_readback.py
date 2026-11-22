"""Independent saved-state checks of the full connection-only radial functional."""
import json,pickle
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
D=pickle.loads((HERE/'normal_ghostjet_insertion.pkl').read_bytes())
R=pickle.loads((HERE/'normal_connection_radial_reduction.pkl').read_bytes())
x=R['radial_variable'];U=R['upper_radial_limit'];ks=D['green_symbols'];hs=D['h_symbols'];P=D['projection']
E=D['on_shell_angular_numerators']['combined']
rho=next(v for v in E.free_symbols if str(v)=='rho');g=next(v for v in E.free_symbols if str(v)=='g');d=x+rho**2
t=s.Symbol('dimensionless_radius',positive=True)
replace={x:rho**2*t,hs[0]:rho**2*hs[0],P:g**2*rho**2*P}
for off,dimension in ((0,0),(3,1),(6,0),(9,0)):
 for j in range(3):replace[ks[off+j]]=g**2*ks[off+j]/rho**(2*(dimension+j))
assert s.expand(E.subs(replace,simultaneous=True)/(g**2*rho**26)-E.subs({rho:1,g:1,x:t},simultaneous=True))==0
literal_h=s.diff(E,hs[0])-g**2*rho**2*(s.diff(E,ks[0])+s.diff(E,ks[6]))/(16*d**2)+g**2*rho**2*(s.diff(E,ks[1])+s.diff(E,ks[7]))/(8*d**3)
homogeneous=s.diff(E,ks[6])/d-s.diff(E,ks[7])/d**2+4*s.diff(E,P)
assert s.factor(literal_h)==0 and s.factor(homogeneous)==0
k,kp,kpp=s.symbols('k kprime ksecond')
infinity_homogeneous=0
for row in R['adjoints']:
 off=row['offset'];l=row['l'];c=row['c'];A=row['A'];B=row['B'];W=row['weight'];z=row['adjoint'];f=row['forcing'];V=4*c/(1+x)-8*x/(1+x)**2
 assert s.cancel(A-s.diff(B,x)-W)==0
 assert s.cancel(4*x*s.diff(z,x,2)-4*l*s.diff(z,x)+V*z-W)==0
 boundary=B*k-4*x*(z*kp-s.diff(z,x)*k)-4*(l+1)*z*k
 derivative=s.diff(boundary,x)+s.diff(boundary,k)*kp+s.diff(boundary,kp)*kpp
 assert s.cancel((derivative+z*f-A*k-B*kp).subs(kpp,(f-4*(l+2)*kp-V*k)/(4*x)))==0
 assert s.limit(B,x,0)==0 and s.limit(z,x,0)==0 and s.limit(x*s.diff(z,x),x,0)==0
 # Saved literal Green sections: k1=O(1/x), k3=O(1/x^2), k_eta1=O(log(x)/x^2).
 profile=1/x if l==1 else (s.log(x)/x**2 if l==2 else 1/x**2)
 assert s.limit(boundary.subs({k:profile,kp:s.diff(profile,x)}),x,s.oo)==0
 if l==2:
  infinity_homogeneous=s.limit(boundary.subs({k:1/(1+x),kp:-1/(1+x)**2}),x,s.oo)
assert s.simplify(s.diff(R['unit_finite_bulk'],U)-R['reduced_unit_integrand'].subs(x,U))==0
assert s.limit(R['unit_finite_bulk'],U,0)==0
compact=s.Poly(R['compact_polynomial'],R['compact_variable'],R['compact_log_symbol'])
independent=sum(c*(-s.Rational(1,(j+1)**2) if m else s.Rational(1,j+1)) for (j,m),c in compact.terms())
assert s.factor(independent-R['unit_whole_integral'])==0
coefficient=s.Rational(9042432,27941760);constant=s.Rational(724901,27941760)
expected=-s.pi**2*(coefficient*P/rho**4+constant*g**2/rho**2)
assert s.factor(R['restored_whole_integral']-expected)==0
finite_P=s.factor(s.diff(R['unit_finite_bulk'],P))
assert s.cancel(s.diff(finite_P,U)-16*U**4*(51*U-463)/(9*(1+U)**12))==0
assert finite_P.subs(U,0)==0 and s.limit(finite_P,U,s.oo)==-coefficient
assert s.limit((finite_P+coefficient)*s.log(U),U,s.oo)==0
log_coefficient=-coefficient/16
finite_part=coefficient*s.Rational(7,48)-constant
# A retained homogeneous k_eta1=C/(1+x) has a nonzero upper face.
# It precisely cancels the shift P -> P+4C; the whole-line literal-section formula must not silently drop it.
assert s.factor(infinity_homogeneous-4*coefficient)==0
print('Full scaling, h0, homogeneous invariance, all adjoints, exact finite faces and independent moments: PASS')
print('Whole-position at fixed P_L, literal Green representative:',expected)
print('Finite-ball log coefficient:',log_coefficient,'finite part:',finite_part,'in units pi² g²/rho²')
print('Homogeneous eta1 upper face:',infinity_homogeneous,'per unit added C')
report={'exact_scaling':True,'literal_h0_invariance':True,'pointwise_allowed_homogeneous_invariance':True,'all_three_adjoint_ODEs':True,'finite_boundary_identity_including_B_k':True,'literal_Green_endpoint_limits':True,'independent_compact_moments':True,'whole_position_literal_section_fixed_P':str(expected),'finite_projection_coefficient':str(finite_P),'finite_ball_log_coefficient':str(log_coefficient),'finite_ball_finite_part':str(finite_part),'homogeneous_eta1_infinity_face_per_C':str(infinity_homogeneous),'physical_scalar_coefficient':'OPEN','scope':'full source-weighted connection-only M1, all b+ghost jets; other physical/BV/composite sectors excluded'}
(HERE/'normal_connection_radial_readback.json').write_text(json.dumps(report,indent=2)+'\n')
