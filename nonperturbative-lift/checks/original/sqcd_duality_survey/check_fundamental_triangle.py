from pathlib import Path
import sympy as s
import json,hashlib
def derivative(term,var):
    coef,odds,evens=term
    odds=list(odds);evens=dict(evens)
    if var.startswith('g'):
        if evens.get(var,0)==0:return (0,(),{})
        coef*=evens[var];evens[var]-=1
    else:
        if var not in odds:return (0,(),{})
        j=odds.index(var);coef*=(-1)**j;odds.pop(j)
    return coef,tuple(odds),evens
def edge(term,outer,inner):
    term=derivative(derivative(term,inner),outer)
    return -term[0],term[1],term[2]
initial=(1,('c1','b1','c2','b2','b3'),{'g1':1,'g2':1,'g3':1})
cycles=[
[('b1','g2'),('g1','b3'),('b2','g3')],
[('g1','b2'),('b1','g3'),('g2','b3')]
]
signs=[]
for cycle in cycles:
    term=initial
    for outer,inner in reversed(cycle):term=edge(term,outer,inner)
    assert term[1]==('c1','c2') and not any(term[2].values())
    signs.append(term[0])
v1,v2,t1,t2,dv1,dv2,dt1,dt2,x,y,z=s.symbols('v1 v2 t1 t2 dv1 dv2 dt1 dt2 x y z')
v=s.Matrix([v1,v2]);t=s.Matrix([[t1,t2]])
dv=s.Matrix([dv1,dv2]);dt=s.Matrix([[dt1,dt2]])
G=s.Matrix([[v1,-t2],[v2,t1]])
dG=s.Matrix([[dv1,-dt2],[dv2,dt1]])
dc=s.Matrix([[x,y],[z,-x]])
h=(t*v)[0]
polytrace=s.trace(dG*G.adjugate()*dc)
polynative=(t*dc*dv-dt*dc*v)[0]
factor=s.Rational(1,2)*(-s.Rational(1,8))*s.Rational(1,2)*s.Rational(1,4)*(-2)
checks={
'cycle_left_derivative_signs':signs==[1,1],
'prefactor_without_g2_pi2':factor==s.Rational(1,64),
'frame_determinant':s.expand(G.det()-h)==0,
'polynomial_hK_formula':s.expand(polytrace-polynative)==0,
'completion_coefficient_sign':s.simplify(-3+s.I*3*(-s.I))==0,
}
assert all(checks.values())
p=Path(__file__).with_name('fundamental_triangle_coefficient_audit.md')
result={'checks':checks,'cycle_signs':signs,'prefactor':str(factor),'note_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
Path(__file__).with_name('fundamental_triangle_readback.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
