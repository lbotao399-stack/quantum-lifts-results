from pathlib import Path
import sympy as s
import json, hashlib
q=s.Matrix(2,4,lambda a,i:s.Symbol(f'q{a}{i}'))
t=s.Matrix(4,2,lambda j,a:s.Symbol(f't{j}{a}'))
M=s.Matrix(4,4,lambda i,j:s.Symbol(f'M{i}{j}'))
mu,m,L=s.symbols('mu m L',nonzero=True)
W=sum(M[i,j]*sum(t[j,a]*q[a,i] for a in range(2)) for i in range(4) for j in range(4))/mu+m*M[3,3]
linear=sum(q[a,i]*s.diff(W,q[a,i]) for i in range(3) for a in range(2))-sum(M[i,j]*s.diff(W,M[i,j]) for i in range(3) for j in range(4))
h=(t[3,:]*q[:,3])[0]
P=s.eye(2)-q[:,3]*t[3,:]/h
pq=P*q[:,:3]
projected=sum(pq[a,i]*s.diff(W,q[a,i]) for i in range(3) for a in range(2))-sum(M[i,j]*s.diff(W,M[i,j]) for i in range(3) for j in range(3))+sum(sum(M[i,j]*(t[3,:]*q[:,i])[0]/h for i in range(3))*s.diff(W,M[3,j]) for j in range(3))
det=M[:3,:3].det()
Winst=-L*det/(mu**4*m)
euler=-sum(M[i,j]*s.diff(Winst,M[i,j]) for i in range(3) for j in range(3))
checks={
'linear_uv_variation_zero':s.expand(linear)==0,
'projector_idempotent':all(s.cancel(x)==0 for x in P*P-P),
'projector_kills_heavy_q':all(s.cancel(x)==0 for x in P*q[:,3]),
'projector_kills_heavy_antiquark':all(s.cancel(x)==0 for x in t[3,:]*P),
'projector_trace_one':s.cancel(s.trace(P)-1)==0,
'rational_uv_variation_zero':s.cancel(projected)==0,
'instanton_source_derivative':s.simplify(L*s.diff(Winst,L)-Winst)==0,
'euler_coefficient_and_sign':s.simplify(euler-3*L*det/(mu**4*m))==0,
'ward_matches_source_insertion':s.simplify(euler+3*L*s.diff(Winst,L))==0,
}
assert all(checks.values()),checks
p=Path(__file__).with_name('microscopic_konishi_route.md')
result={'checks':checks,'note_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
out=Path(__file__).with_name('microscopic_konishi_readback.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
