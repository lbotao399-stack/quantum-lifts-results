"""Exact checks of the source's one-instanton bosonic integral."""
import json
import sympy as s

rho,A,eps,R,t,alpha=s.symbols('rho A eps R t alpha',positive=True)
tail=(1-(eps**2/(eps**2+A))**5)/(10*A)
assert s.simplify(s.diff(tail,eps)+eps**9/(eps**2+A)**6)==0
assert s.limit(tail,eps,s.oo)==0
assert s.simplify(s.limit(tail,eps,0)-1/(10*A))==0
feynman_factor=s.factorial(7)/(s.factorial(3)**2)
x_factor=s.pi**2*s.factorial(5)/s.factorial(7)
assert feynman_factor*x_factor==10*s.pi**2/3
mass=s.pi**2/s.Integer(3)*s.integrate(alpha**2*(1-alpha)**2,(alpha,0,1))
assert mass==s.pi**2/90
assert 36*mass==2*s.pi**2/5
moment=s.Rational(1,2)*(s.log((R**2+A)/A)-sum((R**2/(R**2+A))**j/s.Integer(j) for j in range(1,6)))
assert s.simplify(s.diff(moment,R)-R**11/(R**2+A)**6)==0
assert s.limit(moment,R,0)==0
assert s.limit(s.diff(moment,R)*R,R,s.oo)==1
moment_log_R_coefficient=10*s.pi**2/s.Integer(3)*s.integrate(alpha**3*(1-alpha)**3,(alpha,0,1))
assert moment_log_R_coefficient==s.pi**2/42
f_cubed_integral=s.pi**2*s.factorial(3)/s.factorial(5)/rho**2
assert f_cubed_integral==s.pi**2/(20*rho**2)
print(json.dumps({'status':'PASS','moduli_mass':str(mass),'with_grassmann_factor_36':str(36*mass),
 'rho_cutoff_tail':str(tail),'rho2_moment_log_R_coefficient_divided_by_r2':str(moment_log_R_coefficient),
 'integral_f_cubed':str(f_cubed_integral)},indent=2))
