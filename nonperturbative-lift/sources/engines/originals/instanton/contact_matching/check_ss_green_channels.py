"""Exact angular channels and convergent radial Green construction for SS response.

This verifies the two covariant scalar channels, their regular homogeneous
solutions, the leading Green coefficients, and the match to the independently
constructed free leading response after the singular-gauge rotation.
"""
import sys,pickle
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from physical_profiles import y,rho,den,covariant_gradient
from external_fermion import rotation_numerator as Y
data=pickle.loads((HERE/'logarithmic_connection.pkl').read_bytes())
leading=pickle.loads((HERE/'fermion_asymptotic_break.pkl').read_bytes())
g=next(value for value in data['q'].free_symbols if str(value)=='g')
radial=sum(value**2 for value in y)
ss1=tuple(s.cancel(data['dirac_source']['h'][i][1,j]*den**4/(g**2*rho**4))
          for i in range(3) for j in range(2))
ss3=tuple(s.cancel(data['dirac_source']['hprime'][i][1,j]*den**3/(g**2*rho**2))
          for i in range(3) for j in range(2))
P1=s.Matrix(3,2,ss1)
P3=s.Matrix(3,2,ss3)
H3=(P3-radial*P1/3).applyfunc(s.expand)
P2=s.Matrix(3,1,lambda i,j:s.cancel(data['dirac_source']['h'][i][1,2]*den**4/(g**2*rho**4)))
assert all(data['dirac_source']['hprime'][i][:,2]==s.zeros(2,1) for i in range(3))
assert all(s.expand(sum(value.diff(x,2) for x in y))==0 for value in P1)
assert all(s.expand(sum(value.diff(x,2) for x in y))==0 for value in H3)
assert all(s.expand(sum(value.diff(x,2) for x in y))==0 for value in P2)
print('SS source angular decomposition into harmonic degrees one and three: PASS',flush=True)

# The supplied connection has y.A=0, so radial cross-derivatives are ordinary.
for l,c,poly in ((1,2,P1),(3,5,H3),(2,4,P2)):
    fields=tuple(poly.row(i) for i in range(3))
    laplace=[s.zeros(1,poly.cols) for _ in range(3)]
    for mu in range(4):
        first=covariant_gradient(fields,mu)
        second=covariant_gradient(first,mu)
        for i in range(3):laplace[i]+=second[i]
    for i in range(3):
        assert all(s.cancel(value)==0 for value in laplace[i]-(4*c/den-8*radial/den**2)*fields[i])
    print('Covariant scalar angular channel:',(l,c),'PASS',flush=True)

rad=s.Symbol('s',positive=True)
for l,c,u in ((1,2,(rad+3*rho**2)/(rad+rho**2)),(3,5,1/(rad+rho**2)),(2,4,1/(rad+rho**2))):
    image=4*rad*s.diff(u,rad,2)+4*(l+2)*s.diff(u,rad)+(4*c/(rad+rho**2)-8*rad/(rad+rho**2)**2)*u
    assert s.cancel(image)==0
print('All regular homogeneous radial solutions: PASS',flush=True)

h=data['h'].subs({y[0]:0,y[1]:0,y[2]:0,y[3]:s.sqrt(rad)})
hp=data['hprime'].subs({y[0]:0,y[1]:0,y[2]:0,y[3]:s.sqrt(rad)})
f1=g**2*rho**4*h/(rad+rho**2)**4+g**2*rho**2*rad*hp/(3*(rad+rho**2)**3)
f3=g**2*rho**2*hp/(rad+rho**2)**3
f2=g**2*rho**4*h/(rad+rho**2)**4
assert s.limit(rad**2*f1,rad,s.oo)==g**2*rho**2/24
assert s.limit(rad**3*f3,rad,s.oo)==g**2*rho**2/8
# For v'=s^(-l-2)u^(-2) int_0^s t^(l+1)u(t)f(t)/4 dt,
# l'Hopital gives s^2 v' -> g^2 rho^2/96 and /32 respectively.
for l,u,forcing,expected in ((1,(rad+3*rho**2)/(rad+rho**2),f1,g**2*rho**2/96),
                            (3,1/(rad+rho**2),f3,g**2*rho**2/32)):
    assert s.limit(rad**(l+1)*u*forcing/4,rad,s.oo)==expected
    assert s.limit(rad**(l+1)*u**2/rad**2,rad,s.oo)==1
print('Exact radial Green leading coefficients, independent of h0: PASS',flush=True)
assert s.limit(rad**3*f2,rad,s.oo)==g**2*rho**4/8
assert s.limit(rad**4*f2/(4*(rad+rho**2)),rad,s.oo)==g**2*rho**4/32
print('Eta-one channel: logarithmic numerator primitive; no leading 1/R response: PASS',flush=True)

regular_aux=-g**2*rho**2*(P1/(96*radial)+H3/(32*radial**2))
singular_aux=(Y.T*regular_aux/radial).applyfunc(s.cancel)
for i in range(3):
    for j in range(2):
        assert s.cancel(singular_aux[i,j]-leading['auxiliary_leading'][i][1,j])==0
        assert leading['auxiliary_leading'][i][0,j]==0
print('Global regular Green leading scalar equals independent singular-gauge free inverse: PASS',flush=True)
(HERE/'ss_green_channels.pkl').write_bytes(pickle.dumps({'P1':P1,'P3_harmonic':H3,
    'radial_variable':rad,'forcing_l1':f1,'forcing_l3':f3,
    'P2_eta1':P2,'forcing_l2_eta1':f2,
    'regular_auxiliary_leading':regular_aux,'singular_auxiliary_leading':singular_aux}))
