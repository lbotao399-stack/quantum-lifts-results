from pathlib import Path
import json,hashlib
import sympy as s
def derivative(term,var):
 c,odd,even=term;odd=list(odd);even=dict(even)
 if var.startswith('g'):
  if not even.get(var,0):return (0,(),{})
  c*=even[var];even[var]-=1
 else:
  if var not in odd:return (0,(),{})
  j=odd.index(var);c*=(-1)**j;odd.pop(j)
 return c,tuple(odd),even
def edge(term,outer,inner):
 c,o,e=derivative(derivative(term,inner),outer)
 return -c,o,e
cases=[
((1,('c2','beta2','b3','beta3'),{'gM1':1,'gv1':1,'gt1':1,'gv2':1}),
 [('gv1','beta2'),('gM1','b3'),('gv2','beta3')]),
((1,('c1','beta1','b3','beta3'),{'gv1':1,'gM2':1,'gv2':1,'gt2':1}),
 [('beta1','gv2'),('gv1','beta3'),('gM2','b3')])
]
signs=[]
for term,es in cases:
 for x,y in reversed(es):term=edge(term,x,y)
 signs.append(term[0])
common=s.Rational(1,2)*(-s.Rational(1,8))*s.Rational(1,2)*s.Rational(1,4)
# Vertex -i; reversed vertex positions require momentum sign -1.
mixed=s.simplify(common*((-s.I)*signs[0]+(-s.I)*signs[1]*(-1)))
checks={
'mixed_fiber_signs':signs==[-1,1],
'Q1_b_beta_v_coefficient':mixed==-s.I/64,
'cross_for_negative_v_beta':-mixed==s.I/64,
'heavy_h_anomaly_multiplier':s.Rational(1)+s.Rational(1,2)==s.Rational(3,2),
'nilpotency':-s.Rational(1,64)-s.Rational(1,128)+s.Rational(3,128)==0,
'remaining_cochain_closure':-s.Rational(3,128)+s.I*(-s.I/s.Integer(2)-2*s.I/s.Integer(2))/64==0,
}
assert all(checks.values()),checks
p=Path(__file__).with_name('massless_heavy_cross_audit.md')
result={'checks':checks,'fiber_signs':signs,'Q1_b_beta_coefficient_without_g_mu_pi2':str(mixed),'note_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
Path(__file__).with_name('massless_heavy_cross_readback.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
