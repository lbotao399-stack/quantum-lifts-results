"""Small exact checks of the displayed SQCD F terms; no numerical fitting."""
import json
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
V={(i,j):s.Symbol(f'V{i}{j}') for i in range(1,7) for j in range(i+1,7)}
def pf(indices):
 if not indices:return s.Integer(1)
 return s.expand(sum((-1)**(j+1)*V[indices[0],indices[j]]*pf(indices[1:j]+indices[j+1:]) for j in range(1,len(indices))))
full=pf(tuple(range(1,7)))
for i,j in V:
 assert s.expand(s.diff(full,V[i,j])-(-1)**(i+j+1)*pf(tuple(k for k in range(1,7) if k not in (i,j))))==0
assert s.diff(full,V[1,2])==V[3,4]*V[5,6]-V[3,5]*V[4,6]+V[3,6]*V[4,5]
zero=dict.fromkeys(V.values(),0)
assert all(s.diff(full,v).subs(zero)==0 for v in V.values())
assert all(s.diff(full,v,w).subs(zero)==0 for v in V.values() for w in V.values())
X,M11,M12,M21,M22,B,Bt,m2,L4=s.symbols('X M11 M12 M21 M22 B Btilde m2 Lambda4',nonzero=True)
W=X*(M11*M22-M12*M21-B*Bt-L4)+m2*M22
section={X:-m2/M11,M12:0,M21:0,M22:L4/M11,B:0,Bt:0}
for v in (X,M12,M21,M22,B,Bt):assert s.cancel(s.diff(W,v).subs(section))==0
assert s.cancel(W.subs(section)-m2*L4/M11)==0
assert s.cancel(s.diff(W,M11).subs(section)+m2*L4/M11**2)==0
gamma,L5,m1=s.symbols('gamma Lambda5 m1',nonzero=True)
assert -s.diff(L5/gamma,gamma)==L5/gamma**2
assert s.cancel(gamma**2*(-s.diff(L5/gamma,gamma))-L5)==0
assert s.cancel(-s.diff(L5/gamma+m1*gamma,gamma)-(L5/gamma**2-m1))==0
report={'all_15_six_by_six_Pfaffian_derivatives':True,'displayed_V12_F_term':True,'cubic_origin_gradient_and_Hessian_zero':True,'all_heavy_F_terms_after_exact_elimination':True,'ADS_superpotential_and_remaining_light_F_term':True,'source_sign_Koszul_differential':True,'unit_primitive_after_scale_localization':True,'scope':'small effective SQCD F-term/Koszul calculation; no claim of full electric/magnetic jet-complex equivalence'}
(HERE/'n1_duality_readback.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
