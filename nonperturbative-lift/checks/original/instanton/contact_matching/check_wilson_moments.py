"""Exact positive Wilson-line bound; no floating-point quadrature."""
import sympy as s

t=s.Symbol('t',real=True)
r=s.Symbol('r',positive=True)
normalization=s.integrate(t**2*(1-t**2)**4,(t,0,1))
assert normalization==s.Rational(128,3465)
second=s.integrate(t**4*(1-t**2)**4,(t,0,1))/normalization
assert second==s.Rational(3,13)

def convolution(p,q):
    return s.factor(s.pi**7*r**(14-2*p-2*q)
      *s.gamma(p+q-7)*s.gamma(7-p)*s.gamma(7-q)
      /(s.gamma(p)*s.gamma(q)*s.gamma(14-p-q)))

base=convolution(s.Integer(4),s.Integer(4))
mixed=convolution(s.Rational(7,2),s.Rational(9,2))
both=convolution(s.Rational(9,2),s.Rational(9,2))
assert base==s.pi**7/(1080*r**2)
assert mixed==s.pi**7/(1050*r**2)
assert both==2*s.pi**7/(3675*r**4)
cos_moment=s.factor((2*mixed-r**2*both)/(2*base))
assert cos_moment==s.Rational(36,49)
trace_bound=s.factor(3-2*s.pi**2*second*(1-cos_moment))
assert s.simplify(trace_bound-(3-6*s.pi**2/49))==0
wilson_bound=s.factor(-s.Rational(1,4)+s.Rational(4,27)*trace_bound)
assert s.expand(wilson_bound-(s.Rational(7,36)-8*s.pi**2/441))==0
rational_bound=s.factor(wilson_bound.subs(s.pi,s.Rational(22,7)))
assert rational_bound==s.Rational(1319,86436)>0
assert trace_bound.subs(s.pi,s.Rational(22,7))>s.Rational(27,16)

# Angular variable change yields a fully elementary normalized probability law.
delta=s.Symbol('delta',positive=True)
theta=s.Symbol('theta',real=True)
integrand=s.expand_trig(s.cos(theta)**5*s.cos(theta+delta)**5)
# Direct integration after shifting theta to center its interval reduces to
# powers of cos(delta) and an even angular integral; checked symbolically below.
z=s.Symbol('z',real=True)
half=(s.pi-delta)/2
angular=s.expand((s.cos(delta)+s.cos(2*z))**5/32)
angular_integral=s.integrate(angular,(z,-half,half))
angular_integral=s.trigsimp(s.expand_trig(angular_integral))
density=s.expand_trig(s.sin(delta)*angular_integral)
angular_normalization=s.integrate(density,(delta,0,s.pi))
assert s.simplify(angular_normalization-77*s.pi/1024)==0
angular_cos=s.integrate(s.expand_trig(s.cos(delta)*density),(delta,0,s.pi))
assert s.simplify(angular_cos/angular_normalization-s.Rational(36,49))==0
print('Beta normalization and independence moment: PASS')
print('14-dimensional convolution constants and E[cos Delta]=36/49: PASS')
print('Wilson/scalar >= 7/36 - 8*pi^2/441 > 1319/86436 > 0: PASS')
print('Angular density integral =',angular_normalization)
print('Unnormalized angular density =',s.trigsimp(s.expand_trig(density)))

fourier_density=((s.pi-delta)*(75*s.sin(2*delta)+24*s.sin(4*delta)+s.sin(6*delta))
 +(2000-375*s.cos(2*delta)-1488*s.cos(4*delta)-137*s.cos(6*delta))/60)/(77*s.pi)
assert s.trigsimp(s.expand_trig(fourier_density-density/angular_normalization))==0

def angular_cos_integral(power,k):
    if power==0:return s.Integer(0)
    return -s.Rational(power,k)*angular_sin_integral(power-1,k)
def angular_sin_integral(power,k):
    if power==0:return s.Integer(0)
    return -s.pi**power/k+s.Rational(power,k)*angular_cos_integral(power-1,k)
def angle_moment(power):
    weighted=sum(c*(s.pi*angular_sin_integral(power,k)-angular_sin_integral(power+1,k))
                 for k,c in ((2,75),(4,24),(6,1)))
    elementary=(2000*s.pi**(power+1)/(power+1)
                -sum(c*angular_cos_integral(power,k) for k,c in ((2,375),(4,1488),(6,137))))/60
    return s.expand((weighted+elementary)/(77*s.pi))
assert angle_moment(0)==1
for n in range(1,4):
    direct=s.integrate(s.expand_trig(delta**(2*n)*density),(delta,0,s.pi))/angular_normalization
    assert s.simplify(direct-angle_moment(2*n))==0
    print('E[Delta^%d] ='%(2*n),angle_moment(2*n))
    beta=s.rf(s.Rational(3,2),n)/s.rf(s.Rational(13,2),n)
    assert beta==s.integrate(t**(2*n+2)*(1-t**2)**4,(t,0,1))/normalization
print('Elementary angle moments and rational Beta coefficients for exact series: PASS')
