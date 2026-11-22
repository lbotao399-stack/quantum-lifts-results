import json
from pathlib import Path
import sympy as s
x,u,L,m=s.symbols('x u L m')
curves={
0:x**3-u*x**2+L**4*x/4,
1:x**2*(x-u)-L**6/64,
2:(x**2-L**4/64)*(x-u),
3:x**2*(x-u)-L**2*(x-u)**2/64,
}
expected={
0:L**8*(u**2-L**4)/16,
1:-L**6*(u**3+s.Rational(27,256)*L**6)/16,
2:L**4*(u**2-L**4/64)**2/16,
3:-L**2*u**4*(u-L**2/256)/16,
}
for n in curves:
 assert s.factor(s.discriminant(curves[n],x)-expected[n])==0
# Invert the first two nonconstant terms in a(u) without floating arithmetic.
t=s.symbols('t',positive=True)
c1,c2=s.symbols('c1 c2')
ua=1/(2*t**2)+c1*t**2+c2*t**6
period=s.sqrt(2*ua)*(1-s.Rational(1,16)/ua**2-s.Rational(15,1024)/ua**4)
series=s.series(t*period,t,0,10).removeO()
sol=s.solve([series.coeff(t,4),series.coeff(t,8)],[c1,c2],dict=True)[0]
assert sol=={c1:s.Rational(1,4),c2:s.Rational(5,64)}
# Argyres-Seiberg source coordinate transformation, eqs. (2.3)-(2.4).
xt,yt,vt=s.symbols('xt yt vt',nonzero=True)
expr=(4*yt**2/vt**2-(-(2*s.I*vt)*(2*(-s.I*xt/vt)**3-2*s.I*vt)))*vt**2/4
assert s.simplify(expr-(yt**2-xt**3+vt**4))==0
# Nf4 with masses (m,m,0,0), SW (16.26): coupling-coordinate invariant degeneration data.
e1,e2=s.symbols('e1 e2')
e=[e1,e2,-e1-e2]
p4=s.prod(x-z*u-z*z*m*m for z in e)
d4=s.prod((e[i]-e[j])**2 for i in range(3) for j in range(i+1,3))*s.prod((u-z*m*m)**2 for z in e)
assert s.factor(s.discriminant(p4,x)-d4)==0
out={
 'massless_discriminants':{str(n):str(s.factor(expected[n])) for n in curves},
 'period_coefficients':[str(s.binomial(s.Rational(1,2),2*k)*s.factorial(2*k)/(2**(2*k)*s.factorial(k)**2)) for k in range(4)],
 'inverse_period_coefficients':{str(k):str(v) for k,v in sol.items()},
 'prepotential_F1':str(s.Rational(1,4)/(2*s.pi*s.I)),
 'prepotential_F2':str(s.Rational(5,64)/(4*s.pi*s.I)),
 'Nf3_physical_singularities_in_SW_scale_squared':['-1/432',str(s.Rational(1,256)-s.Rational(1,432))],
 'Nf4_PV_two_instanton_c2':str(-s.Rational(7,2**6*3**5)),
 'Argyres_Seiberg_transform':'passed',
 'Nf4_massive_discriminant':'passed',
 'E6_beta_function':-2*4+2*1+6,
 'Higgs_complex_dimensions':[2*3*6-2*8,2*2+22-2*3]
}
p=Path('work/sqcd_duality_survey/n2_readback.json')
p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
