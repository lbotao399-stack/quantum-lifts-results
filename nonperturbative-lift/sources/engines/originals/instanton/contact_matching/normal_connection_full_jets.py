"""Complete source-weighted ghost-jet connection block, with actual normal fields."""
import sys,pickle,math,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from physical_profiles import y,rho,den,sigma,X,eta,covariant_gradient
from external_fermion import rotation_numerator as Y
N=pickle.loads((HERE/'normal_connection_response.pkl').read_bytes());K=pickle.loads((HERE/'normal_connection_ghost_kernels.pkl').read_bytes());B=pickle.loads((HERE/'normal_bjet_insertion.pkl').read_bytes())
rad=N['radial_variable'];hs=N['h_symbols'];ks=N['green_symbols'];P=B['finite_cutoff_projection'];g=next(v for v in B['angular_average_numerator'].free_symbols if str(v)=='g')
RR,*gen=ring((*y,rho,rad,g,P,*hs,*ks),QQ_I)
yy=gen[:4];rr,ss,gg,pp=gen[4:8];hh=gen[8:13];kk=gen[13:];dd=ss+rr**2;imag=QQ_I.convert(s.I)
def poly(value):return RR.from_expr(s.expand(value))
def partial(value,mu):
 out=value.diff(yy[mu])+2*yy[mu]*value.diff(ss)
 for j in (1,2,3):out+=2*yy[mu]*hh[j+1]*value.diff(hh[j])
 for start in (0,3,6,9):
  for j in (0,1):out+=2*yy[mu]*kk[start+j+1]*value.diff(kk[start+j])
 return out
A=[[poly(2*sum(eta[i][mu][nu]*y[nu] for nu in range(4))) for mu in range(4)] for i in range(3)]
def holo(fields,power,alpha):
 return tuple(sum((QQ_I.convert(sigma[mu][alpha,0])*(dd*partial(fields[i],mu)-2*power*yy[mu]*fields[i]+sum((int(s.LeviCivita(i,j,k))*A[j][mu]*fields[k] for j in range(3) for k in range(3)),RR.zero)) for mu in range(4)),RR.zero) for i in range(3))
H=Y/den;DH=[covariant_gradient(tuple(H.row(i) for i in range(3)),mu) for mu in range(4)]
normal=[[[poly(v) for v in row] for row in color] for color in N['normal_numerators']]
for i in range(3):
 for mu in range(4):normal[i][mu][2]-=pp*dd*poly(s.factor((DH[mu][i][0,0]-s.I*DH[mu][i][0,1])*den**2))
null=(-poly(X[0,0]),-poly(X[1,0]),RR.one,RR.zero)
full_cart={}
for alpha,jet in ((0,(1,0)),(1,(0,1))):
 full_cart[jet]=tuple(tuple(sum((QQ_I.convert(sigma[mu][alpha,0])*normal[i][mu][m] for mu in range(4)),RR.zero) for m in range(4)) for i in range(3))
for jet in ((2,0),(1,1),(0,2)):
 r,t=jet;rows=[[RR.zero for m in range(4)] for i in range(3)]
 for alpha,count,prev in ((0,r,(r-1,t)),(1,t,(r,t-1))):
  if not count:continue
  for m in range(4):
   result=holo(tuple(full_cart[prev][i][m] for i in range(3)),3,alpha)
   for i in range(3):rows[i][m]+=QQ_I.convert(s.Rational(count,2))*result[i]
 full_cart[jet]=tuple(tuple(row) for row in rows)
full_source={}
for jet,row in full_cart.items():
 full_source['E',*jet]=tuple((-imag*row[0][m]-row[1][m])/2 for m in range(4))
 full_source['H',*jet]=tuple(-imag*row[2][m]/2 for m in range(4))
 full_source['F',*jet]=tuple((-imag*row[0][m]+row[1][m])/2 for m in range(4))
old=pickle.loads((HERE/'normal_ghostjet_insertion.pkl').read_bytes())
for key,row in full_source.items():
 assert sum((row[m]*null[m] for m in range(4)),RR.zero)==poly(old['source_normal_jets'][key])
print('All 60 full source normal-jet coefficients agree with their independently differentiated annihilator contraction: PASS',flush=True)
(HERE/'normal_connection_full_jets.pkl').write_bytes(pickle.dumps({'source_normal_jets':{key:tuple(v.as_expr() for v in row) for key,row in full_source.items()},'denominator_power':'total source degree + 2','radial_variable':rad,'h_symbols':hs,'green_symbols':ks,'projection':P}))
