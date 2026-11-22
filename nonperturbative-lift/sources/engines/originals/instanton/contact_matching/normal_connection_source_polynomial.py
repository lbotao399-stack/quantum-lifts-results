"""Source-exact connection variation, preserving full eq-B-14 weights."""
import json,pickle,sys
from pathlib import Path
from fractions import Fraction as F
from math import comb
import sympy as s
HERE=Path(__file__).resolve().parent
raw=json.loads((HERE.parent.parent/'rational_c4/C2_complete_exact.json').read_text())['polynomials']['partner_Nb1']
colors=('E','H','F');basis={'E':s.Matrix([[0,1],[0,0]]),'H':s.diag(1,-1),'F':s.Matrix([[0,0],[1,0]])}
f={}
for left in colors:
 for right in colors:
  bracket=basis[left]*basis[right]-basis[right]*basis[left]
  for output,value in zip(('E','H','F'),(bracket[0,1],bracket[0,0],bracket[1,0])):
   if value:f[output,left,right]=int(value)
def canon(word):
 even=sorted(x for x in word if x[0]!='c');odd=[x for x in word if x[0]=='c']
 if len(set(odd))!=len(odd):return None,0
 sign=(-1)**sum(odd[i]>odd[j] for i in range(len(odd)) for j in range(i+1,len(odd)))
 return tuple(even+sorted(odd)),sign
def image(letter):
 kind,color,r,t=letter;N=r+t;out=[]
 if kind=='b':
  if N==1:
   for zcol in colors:
    for bcol in colors:
     coefficient=-f.get((bcol,zcol,color),0)
     if coefficient:out.append((F(coefficient),(('z',zcol,r,t),('b',bcol,0,0))))
 elif kind=='c' and N>=2:
  for a in range(r+1):
   for b in range(t+1):
    A=a+b
    if not (1<=A<N):continue
    weight=F(N-A,N)*comb(r,a)*comb(t,b)
    for zcol in colors:
     for ccol in colors:
      coefficient=f.get((color,zcol,ccol),0)
      if coefficient:out.append((weight*coefficient,(('z',zcol,a,b),('c',ccol,r-a,t-b))))
 return out
pieces={'bjets':{},'ghostjets':{}}
for term in raw:
 word=tuple((x['kind'],x['color'],*x['jet']) for x in term['letters']);coefficient=F(*term['coefficient'])
 for j,letter in enumerate(word):
  for weight,replacement in image(letter):
   key,sign=canon(word[:j]+replacement+word[j+1:])
   if sign:
    which='bjets' if letter[0]=='b' else 'ghostjets'
    pieces[which][key]=pieces[which].get(key,F(0))+coefficient*weight*sign
pieces={name:{word:c for word,c in poly.items() if c} for name,poly in pieces.items()}
print('Exact source variation term counts:',{name:len(poly) for name,poly in pieces.items()},flush=True)
# Explicit low-jet weights; no full-polarization shortcut is used.
assert image(('c','E',2,0))==[(F(-2),(('z','E',1,0),('c','H',1,0))),(F(2),(('z','H',1,0),('c','E',1,0)))]
# Re-read the original source tree differential over Q.
sys.path.insert(0,str(HERE.parent.parent/'rational_c4'))
import exact_product as source
labels={'E':'e12','H':'h1','F':'e21'}
polynomial={}
for term in raw:
 word=tuple((x['kind']+labels[x['color']],*x['jet']) for x in term['letters'])
 key,sign=source.T.canon(word)
 if sign:polynomial[key]=polynomial.get(key,F(0))+sign*F(*term['coefficient'])
assert source.Q0(polynomial)=={}
q20=source.S.IMG(('ce12',2,0));q30=source.S.IMG(('ce12',3,0));qb10=source.S.IMG(('be12',1,0))
print('Original source Q0 M1=0 over Q: PASS; c20 connection variation has one slot versus two in complete polarization.',flush=True)
(HERE/'normal_connection_source_polynomial.pkl').write_bytes(pickle.dumps({'pieces':pieces,'structure_constants':f,'low_jet_connection_images':{('c','E',2,0):image(('c','E',2,0)),('c','E',3,0):image(('c','E',3,0)),('c','E',2,1):image(('c','E',2,1))},'source_Q0_c20':q20,'source_Q0_c30':q30,'source_Q0_b10':qb10,'source_Q0_M1_zero':True,'scope':'even connection variation after contraction with the verified annihilator; z_(r,t) denotes full symmetric covariant jet of R_alpha contracted with that annihilator'}))
