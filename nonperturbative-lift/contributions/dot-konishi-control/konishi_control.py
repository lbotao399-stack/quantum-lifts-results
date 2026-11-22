"""Independent exterior-algebra reconstruction of AK 1605.06346v2 (2.4)-(2.7).
Only public profile formulas are inputs; NOT a BV-Ward operator computation.
"""
from fractions import Fraction as F
from itertools import permutations,product
from pathlib import Path
import json,time,resource
import sympy as s
start=time.monotonic()
def add(*ps):
 r={}
 for p in ps:
  for m,c in p.items():r[m]=r.get(m,0)+c
 return {m:c for m,c in r.items() if c!=0}
def scale(p,c):return {m:c*v for m,v in p.items() if c*v!=0}
def mul(a,b):
 r={}
 for m,c in a.items():
  for n,d in b.items():
   if m&n:continue
   inv=sum((n&((1<<i)-1)).bit_count() for i in range(m.bit_length()) if m>>i&1)
   k=m|n;r[k]=r.get(k,0)+c*d*(-1)**inv
 return {m:c for m,c in r.items() if c!=0}
def eps(*a):
 a=tuple(map(int,a))
 if len(set(a))!=len(a):return 0
 return (-1)**sum(a[i]>a[j] for i in range(len(a)) for j in range(i+1,len(a)))
def z(A,a):return {1<<(2*A+a):F(1)}
def zero():return [[{} for _ in range(2)] for _ in range(2)]
def madd(a,b):return [[add(a[i][j],b[i][j]) for j in range(2)] for i in range(2)]
def mscale(a,c):return [[scale(a[i][j],c) for j in range(2)] for i in range(2)]
def mmul(a,b):return [[add(*(mul(a[i][k],b[k][j]) for k in range(2))) for j in range(2)] for i in range(2)]
def ptr(a,b):return add(*(mul(a[i][j],b[j][i]) for i,j in product(range(2),repeat=2)))
# Remove common f=8 rho^2/(rho^2+y^2)^2. Fundamental color matrix of Eq 2.7.
Fm={(a,b):s.Matrix(2,2,lambda i,j:eps(i,a)*int(b==j)+eps(i,b)*int(a==j)) for a,b in product(range(2),repeat=2)}
Z2={(A,B):add(*(scale(mul(z(A,a),z(B,b)),eps(a,b)) for a,b in product(range(2),repeat=2))) for A,B in product(range(4),repeat=2)}
ZF={};ZFF={}
for A,B in product(range(4),repeat=2):
 v=zero();w=zero()
 for a,b in product(range(2),repeat=2):
  zz=mul(z(A,a),z(B,b));v=madd(v,[[scale(zz,int(Fm[a,b][i,j])) for j in range(2)] for i in range(2)])
  H=sum((eps(c,d)*Fm[a,c]*Fm[d,b] for c,d in product(range(2),repeat=2)),s.zeros(2))
  w=madd(w,[[scale(zz,int(H[i,j])) for j in range(2)] for i in range(2)])
 ZF[A,B]=v;ZFF[A,B]=w
# phi2=ZF/sqrt(2); phi6= -H6/(20 sqrt(2)); lower phi = epsilon phi /2.
H6={}
for A,B in product(range(4),repeat=2):
 H=zero()
 for C,D,E,G in permutations(range(4)):
  t=mul(Z2[A,C],Z2[B,D]);H=madd(H,[[scale(mul(t,ZFF[E,G][i][j]),eps(C,D,E,G)) for j in range(2)] for i in range(2)])
 H6[A,B]=H
K4={};K8={}
for A,B,C,D in permutations(range(4)):
 K4=add(K4,scale(ptr(ZF[A,B],ZF[C,D]),F(eps(A,B,C,D),4)))
 K8=add(K8,scale(ptr(ZF[A,B],H6[C,D]),F(-eps(A,B,C,D),40)))
assert K4=={},K4
assert K8=={255:F(-576)},K8
assert F(-576)*8**3==-9*2**15
# Independent one-flavor Berezin saturation, direct exterior multiplication vs 4x4 determinant.
a,b,c,d,e,f,g,h=s.symbols('a b c d e f g h')
X=s.Matrix([[a,b],[c,d]]);Y=s.Matrix([[e,f],[g,h]])
forms=[{1:1,4:a,8:b},{2:1,4:c,8:d},{1:1,4:e,8:f},{2:1,4:g,8:h}]
w={0:1}
for v in forms:w=mul(w,v)
wed=s.expand(w[15]);det=s.det(s.BlockMatrix([[s.eye(2),X],[s.eye(2),Y]]).as_explicit())
assert s.expand(wed-det)==0 and s.expand(wed-(Y-X).det())==0
# All 16 modes at fixed noncoincident points: X=0,Y=[[2,1],[-1,3]], det=7.
full={0:1}
for A in range(4):
 for row in [[1,0,0,0],[0,1,0,0],[1,0,2,1],[0,1,-1,3]]:
  full=mul(full,{1<<(4*A+i):v for i,v in enumerate(row) if v})
assert full=={65535:7**4},full
# Grouped-by-flavor measure to all xi followed by all eta is even (24 crossings).
assert sum(2*2*(3-A) for A in range(4))==24
# Source measure g^8/(2^34 pi^10), two normalized K insertions g^-4.
coefficient=F((-9*2**15)**2,2**34)
assert coefficient==F(81,16)
# Independent Schwinger-parameter evaluation of remaining bosonic integral.
# d=4-2eps; t=rho^2. rho^7 drho=t^3 dt/2.
# Gaussian x0 integral, t integral, then beta integral yield:
# I=pi^(2-eps)*Gamma(4)*Gamma(6+eps)/(2 Gamma(6)^2)
#   *Gamma(-eps)^2/Gamma(-2eps) * (x12^2)^(-2-eps).
epsilon=s.symbols('epsilon',real=True)
radial=s.gamma(4)*s.gamma(6+epsilon)*s.gamma(-epsilon)**2/(2*s.gamma(6)**2*s.gamma(-2*epsilon))
residue=F(s.limit(epsilon*radial,epsilon,0)) # coefficient of pi^2/eps, gamma ratio residue -2 cancels /2
assert residue==F(-1,20)
# Born=36/(16 pi^4 x^4). The pole ratio multiplies g^4/pi^4.
ratio=coefficient*residue/F(36,16)
assert ratio==F(-9,80) # equivalently -9/5*(g^2/(4pi^2))^2
out={'scope':'public N=4 Konishi normalization control, not j24 Ward image','K4_polynomial':{},'K8_over_f_cubed':str(K8[255]),'K8_over_rho6_denominator6':-9*2**15,'single_flavor_top':str(s.factor(wed)),'sixteen_mode_test':7**4,'measure_reordering_crossings':24,'two_point_prefactor_g4_pi_minus10':str(coefficient),'bosonic_pole_pi2':str(residue),'normalized_pole_g4_pi_minus4':str(ratio),'source_equations':['2.4','2.5','2.6','2.7','2.10','1.3'],'comparison_equations':['2.16','2.17','2.23','2.24','2.26','2.27'],'profile_input_not_rederived':'phi2 and phi6 are public source solutions, not independently solved field equations','runtime_s':time.monotonic()-start,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_name('konishi-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
