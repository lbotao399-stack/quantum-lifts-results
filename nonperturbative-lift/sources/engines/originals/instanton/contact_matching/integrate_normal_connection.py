"""Exact finite-face adjoint reduction of the full connection-only M1 block."""
import pickle
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
D=pickle.loads((HERE/'normal_ghostjet_insertion.pkl').read_bytes())
G=pickle.loads((HERE/'ss_green_channels.pkl').read_bytes())
x=D['radial_variable'];ks=D['green_symbols'];hs=D['h_symbols'];P=D['projection']
angular=D['on_shell_angular_numerators']['combined']
rho=next(v for v in angular.free_symbols if str(v)=='rho')
g=next(v for v in angular.free_symbols if str(v)=='g')
unit=s.factor_terms(x*angular.subs({rho:1,g:1})/(1+x)**16)
explicit=s.factor(unit.subs(dict.fromkeys(ks,0)))
answer=explicit;records=[]
for off,l,c,forcing in ((0,1,2,G['forcing_l1']),(3,3,5,G['forcing_l3']),(6,2,4,G['forcing_l2_eta1'])):
 A=s.factor(s.diff(unit,ks[off]));B=s.factor(s.diff(unit,ks[off+1]))
 W=s.factor(A-s.diff(B,x))
 coeff=s.symbols('adjoint_1:14')
 trial=sum(v*x**j for j,v in enumerate(coeff,1))/(1+x)**10
 residual=s.cancel(4*x*s.diff(trial,x,2)-4*l*s.diff(trial,x)+(4*c/(1+x)-8*x/(1+x)**2)*trial-W)
 solution=s.linsolve(s.Poly(residual.as_numer_denom()[0],x).all_coeffs(),coeff)
 assert len(solution)==1,(l,solution)
 values=next(iter(solution));free=set().union(*(v.free_symbols for v in values))
 values=tuple(v.subs(dict.fromkeys(free,0)) for v in values)
 z=s.factor(trial.subs(dict(zip(coeff,values))))
 assert s.cancel(4*x*s.diff(z,x,2)-4*l*s.diff(z,x)+(4*c/(1+x)-8*x/(1+x)**2)*z-W)==0
 force=forcing.subs({rho:1,g:1}).subs(G['radial_variable'],x)
 answer+=z*force
 records.append({'offset':off,'l':l,'c':c,'A':A,'B':B,'weight':W,'adjoint':z,'forcing':force})
 print('Adjoint',l,':',z,flush=True)
log=s.Symbol('logarithm')
answer=s.cancel(answer.subs(s.log(1+x),log))
num,den=answer.as_numer_denom()
power=next(n for f,n in s.factor_list(den,x)[1] if s.expand(f-(1+x))==0)
polynomial=s.Poly(s.cancel(answer*(1+x)**power),x,log)
integral=0
for (j,k),coefficient in polynomial.terms():
 assert 0<=j<=power-2 and k<=1,(j,k,power)
 beta=s.factorial(j)*s.factorial(power-j-2)/s.factorial(power-1)
 integral+=coefficient*beta*(s.harmonic(power-1)-s.harmonic(power-j-2) if k else 1)
integral=s.factor(integral)
answer=answer.subs(log,s.log(1+x))
u=s.Symbol('u',positive=True);U=s.Symbol('upper_radial_limit',positive=True)
compact=s.Poly(s.cancel(answer.subs(s.log(1+x),-s.log(u)).subs(x,(1-u)/u)/u**2).expand().subs(s.log(u),log),u,log)
independent=0;finite=0
for (j,k),coefficient in compact.terms():
 assert j>=0 and k<=1,(j,k)
 independent+=coefficient*(-s.Rational(1,(j+1)**2) if k else s.Rational(1,j+1))
 lower=1/(1+U)
 primitive=lower**(j+1)*(s.log(lower)/(j+1)-s.Rational(1,(j+1)**2)) if k else lower**(j+1)/(j+1)
 upper=-s.Rational(1,(j+1)**2) if k else s.Rational(1,j+1)
 finite+=coefficient*(upper-primitive)
assert s.factor(independent-integral)==0
assert s.simplify(s.diff(finite,U)-answer.subs(x,U))==0
assert s.limit(finite,U,0)==0
restored=s.factor(s.pi**2*g**2/rho**2*integral.subs({hs[0]:hs[0]/rho**2,P:P/(g**2*rho**2)}))
print('Whole-position at fixed P_L:',restored,flush=True)
saved={'radial_variable':x,'unit_integrand':unit,'adjoints':records,'reduced_unit_integrand':answer,'compact_polynomial':compact.as_expr(),'compact_variable':u,'compact_log_symbol':log,'unit_whole_integral':integral,'restored_whole_integral':restored,'unit_finite_bulk':finite,'upper_radial_limit':U,'finite_boundary_formula':'sum_l [B_l k_l -4U(z_l k_lprime-z_lprime k_l)-4(l+1)z_l k_l] at U; finite integral equals bulk plus this sum','scope':'full source-weighted connection-only M1, b+ghost jets; finite P_L and all boundary forms retained; physical scalar coefficient OPEN'}
(HERE/'normal_connection_radial_reduction.pkl').write_bytes(pickle.dumps(saved))
