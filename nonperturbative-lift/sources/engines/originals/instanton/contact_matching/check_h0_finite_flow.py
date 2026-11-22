"""Exact finite coordinate flow of the literal h0 section vector.

This integrates the verified linear-response section vector; it does not claim
that the higher-background physical field family is known.
"""
import pickle,json
from pathlib import Path
import sympy as s
def factor(value):return s.cancel(value)
HERE=Path(__file__).resolve().parent
D=pickle.loads((HERE/'h0_color_measure.pkl').read_bytes())
F=pickle.loads((HERE/'h0_color_fermions.pkl').read_bytes())
n=D['n'];rho=next(v for v in D['scale_coefficient'].free_symbols if str(v)=='rho')
g=next(v for v in D['scale_coefficient'].free_symbols if str(v)=='g')
t,u=s.symbols('section_parameter other_section_parameter')
z=n[0]-s.I*n[1];a=g**2*z;den=1-t*a
omega=D['orientation_coefficients']
cross=s.Matrix(3,3,lambda i,j:sum(s.LeviCivita(i,k,j)*omega[k] for k in range(3)))
rotation=s.eye(3)+t*(2-t*a)*cross/(2*den)+t**2*cross**2/(2*den)
rotation=rotation.applyfunc(factor)
newn=(rotation*s.Matrix(n)).applyfunc(factor)
assert all(s.simplify(factor(v))==0 for v in rotation.T*rotation-s.eye(3))
assert factor(rotation.det()-1)==0
assert factor((newn.T*newn)[0]-sum(v*v for v in n))==0
assert factor(newn[0]-s.I*newn[1]-z/den)==0
assert factor(newn[2]-n[2]/den)==0
for i in range(3):
 assert factor(s.diff(newn[i],t)-sum(s.LeviCivita(i,j,k)*omega[j].subs(dict(zip(n,newn)),simultaneous=True)*newn[k] for j in range(3) for k in range(3)))==0
odd=s.eye(4);odd[2,2]=1/den;odd[3,3]=den;odd[3,2]=-g**2*n[2]*t*(2-t*a)/(2*den)
C=F['fermion_coordinate_matrix']
assert all(s.simplify(factor(v))==0 for v in odd.diff(t)-C.subs(dict(zip(n,newn)),simultaneous=True)*odd)
assert factor(odd.det()-1)==0
# rho maps to rho*sqrt(den), on the branch equal to one at t=0.
assert factor(s.diff(rho*s.sqrt(den),t)-D['scale_coefficient'].subs(dict(zip(n,newn))|{rho:rho*s.sqrt(den)},simultaneous=True))==0
print('Exact framing/color, size and four-fermion flow equations: PASS',flush=True)
# Complete three-dimensional framing differential from its variational ODE.
null=s.Matrix([1,-s.I,0]);norm=sum(v*v for v in n)
frame_jacobian=(s.eye(3)-g**2*t*null*s.Matrix(n).T+g**4*t**2*norm*null*null.T/2)/den
variational=(a/den)*s.eye(3)-g**2*null*newn.T
assert all(factor(v)==0 for v in frame_jacobian.diff(t)-variational*frame_jacobian)
assert frame_jacobian.subs(t,0)==s.eye(3)
haar=factor(frame_jacobian.det())
assert factor(haar-den**-2)==0,haar
plain_super=s.sqrt(den)*haar/odd.det()
weighted=factor(den*haar/odd.det())
assert factor(weighted-1/den)==0
assert factor(s.diff(s.log(weighted),t)-a/den)==0
print('Finite Haar determinant = (1-t g² z)^(-2); odd determinant = 1; rho-weighted Berezinian = (1-t g² z)^(-1): PASS',flush=True)
# Verify finite group law for framing and all odd coordinates.
newn_u=newn.subs(t,u)
rot_t_after_u=rotation.subs(t,t/(1-u*a))
rot_u=rotation.subs(t,u)
assert all(s.simplify(factor(v))==0 for v in rot_t_after_u*rot_u-rotation.subs(t,t+u))
odd_after_u=s.eye(4);odd_after_u[2,2]=1/(1-t*a/(1-u*a));odd_after_u[3,3]=1-t*a/(1-u*a);odd_after_u[3,2]=-g**2*n[2]*t*(2-t*a/(1-u*a))/(2*(1-u*a)*(1-t*a/(1-u*a)))
assert all(s.simplify(factor(v))==0 for v in odd_after_u*odd.subs(t,u)-odd.subs(t,t+u))
newden=1-t*a/(1-u*a)
assert factor(den.subs(t,u)*newden-den.subs(t,t+u))==0
assert factor(weighted.subs(t,u)*1/newden-weighted.subs(t,t+u))==0
print('Exact composition law for framing, rho², odd map and density cocycle: PASS',flush=True)
saved={'parameter':t,'other_parameter':u,'n':n,'rho':rho,'g':g,'denominator':den,'framing_rotation':rotation,'finite_n':newn,'finite_rho':rho*s.sqrt(den),'odd_matrix':odd,'framing_jacobian':frame_jacobian,'Haar_jacobian':haar,'plain_coordinate_Berezinian':plain_super,'rho_weighted_Berezinian':weighted,'scope':'exact flow of verified literal h0 coordinate section; finite physical source-family matching beyond linear background order remains OPEN'}
(HERE/'h0_finite_flow.pkl').write_bytes(pickle.dumps(saved))
(HERE/'h0_finite_flow_readback.json').write_text(json.dumps({'flow_equations':True,'framing_Haar_jacobian':'(1-t*g**2*(n1-I*n2))**(-2)','odd_determinant':1,'weighted_Berezinian':'(1-t*g**2*(n1-I*n2))**(-1)','exact_composition':True,'physical_finite_family':'OPEN'},indent=2)+'\n')
