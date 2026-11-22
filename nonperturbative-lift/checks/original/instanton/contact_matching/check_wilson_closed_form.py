"""Exact finite Si/Ci evaluation of the one-instanton adjoint Wilson trace."""
import sympy as s
t=s.Symbol('t',real=True)
delta=s.Symbol('delta',real=True)
for k in (2,4,6):
    for value in (s.Rational(1,3),s.Rational(2,3)):
        direct=s.integrate((s.pi-delta)*s.sin(k*delta)*s.cos(2*value*delta),(delta,0,s.pi))
        target=s.pi*k/(k*k-4*value*value)+4*k*value*s.sin(2*s.pi*value)/(k*k-4*value*value)**2
        assert s.simplify(direct-target)==0
        direct=s.integrate(s.cos(k*delta)*s.cos(2*value*delta),(delta,0,s.pi))
        target=2*value*s.sin(2*s.pi*value)/(4*value*value-k*k)
        assert s.simplify(direct-target)==0
weight=t**2*(1-t**2)**4
R=s.factor(weight*sum(s.Rational(c*k,1)/(k*k-4*t*t)
                     for k,c in ((2,75),(4,24),(6,1))))
S=s.factor(weight*(sum(4*c*k*t/(k*k-4*t*t)**2
                     for k,c in ((2,75),(4,24),(6,1)))
   +s.Rational(50,3)/t
   +sum(s.Rational(c,30)*t/(k*k-4*t*t)
        for k,c in ((2,375),(4,1488),(6,137)))))
R_apart=(-63*t**8+105*t**6-333*t**4-1125*t**2-8088
         +9216/(t+3)+1944/(t+2)-1944/(t-2)-9216/(t-3))
S_apart=(-s.Rational(20736,5)/(t+3)-4608/(t+3)**2
         +s.Rational(20736,5)/(t+2)-972/(t+2)**2
         +s.Rational(20736,5)/(t-2)+972/(t-2)**2
         -s.Rational(20736,5)/(t-3)+4608/(t-3)**2)
assert s.cancel(R-R_apart)==0
assert s.cancel(S-S_apart)==0

R_integral=-s.Rational(42608,5)+9216*s.log(2)+1944*s.log(3)
R_direct=s.integrate(R_apart,(t,0,1))
assert s.simplify(s.expand_log(R_direct,force=True)-R_integral)==0

def single(pole):
    return s.Si(2*s.pi*(1-pole))-s.Si(-2*s.pi*pole)
def double(pole):
    return 2*s.pi*(s.Ci(2*s.pi*abs(1-pole))-s.Ci(2*s.pi*abs(pole)))
S_integral=(s.Rational(20736,5)*(single(-2)+single(2)-single(-3)-single(3))
            +972*(double(2)-double(-2))+4608*(double(3)-double(-3)))
expected=(6480*(s.log(2)+s.Ci(4*s.pi)-s.Ci(8*s.pi))
    +s.Rational(10935,8)*(s.log(3)+s.Ci(2*s.pi)-s.Ci(6*s.pi))
    +2916/s.pi*(s.Si(4*s.pi)+s.Si(6*s.pi)-s.Si(2*s.pi)-s.Si(8*s.pi))
    -s.Rational(23963,4))
assert s.expand(1+s.Rational(45,64)*(R_integral+S_integral/s.pi)-expected)==0

# Verify the definite integral formulas by differentiation at every required
# pole, then use sine(0)=sine(2*pi)=0 for the displayed endpoints.
for pole in (-3,-2,2,3):
    signed=1 if pole<0 else -1
    primitive_single=s.Si(2*s.pi*(t-pole))
    primitive_double=-s.sin(2*s.pi*t)/(t-pole)+2*s.pi*s.Ci(2*s.pi*signed*(t-pole))
    assert s.trigsimp(s.diff(primitive_single,t)-s.sin(2*s.pi*t)/(t-pole))==0
    assert s.trigsimp(s.diff(primitive_double,t)-s.sin(2*s.pi*t)/(t-pole)**2)==0

second_angle=100*s.pi**2/693-s.Rational(91,108)
strong_trace=s.expand(3-4*s.Rational(3,13)*second_angle)
assert strong_trace==s.Rational(34,9)-400*s.pi**2/3003
strong_wilson=s.expand(-s.Rational(1,4)+s.Rational(4,27)*strong_trace)
assert strong_wilson==s.Rational(301,972)-1600*s.pi**2/81081
rational_bound=strong_wilson.subs(s.pi,s.Rational(22,7))
assert rational_bound==s.Rational(5470949,47675628)>0
print('Angular characteristic function rational partial fractions: PASS')
print('Si/Ci primitive differentiation and exact finite expression: PASS')
print('E tr U =',expected)
print('Stronger Wilson/scalar bound > 5470949/47675628 > 0: PASS')
