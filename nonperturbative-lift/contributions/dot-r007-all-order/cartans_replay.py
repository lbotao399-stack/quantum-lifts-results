"""Independent trace-prefactor reconstruction, rational stacking and original minors."""
import json
from pathlib import Path
import sys
from os import environ
from hashlib import sha256
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from functools import lru_cache
import cartan_multiplication as m
vend = Path(environ['DMK_VENDOR']).resolve()
assert sha256((vend/'fortuity/bmn/dmk.py').read_bytes()).hexdigest() == '535ec58677878ffe5b613ec93bc4539c7e58eeefc8e532b98e791ccff0fc86b8'
assert sha256((vend/'fortuity/bmn/words.py').read_bytes()).hexdigest() == 'f74fba895060d3d5fb88922fff4ff38ff8d03546a8d56e7b25700fa19e93fe3e'
sys.path.insert(0,str(vend))
from fortuity.bmn import dmk

@lru_cache(None)
def ev_poly(name, idx):
    def letter(x, row):
        mask=1<<x if x<3 else 7^(1<<(x-3))
        return m.add((((-1)**(x-3) if x>=3 else 1),m.at(mask,row)))
    z=defaultdict(F)
    for multi,c in dmk.tensor(name)(*idx).items():
        P={():F(c)}
        for trace in multi:
            T=[]
            for i in range(3):
                V={():F(1)}
                for x in trace:V=m.mul(V,letter(x,i))
                T.append((1,V))
            P=m.mul(P,m.add(*T))
        for w,a in P.items():z[w]+=a
    return {w:c for w,c in z.items() if c}
got=[{} for i in range(3)]
for c,name,fn in dmk.TERMS:
    if name not in ('s22','s23'):continue
    variables=sorted({a for key,s in fn for a in s.split()}-set('ijk'))
    for inds in product(range(3),repeat=len(variables)):
        I=dict(zip(variables,inds))|dict(i=0,j=0,k=0)
        w=c;P={():F(1)};r=None
        for f,st in fn:
            t=tuple(I[i] for i in st.split())
            if f=='E':w*=dmk.eps(*t)
            elif f=='r16_01':r=t[0]
            else:P=m.mul(P,ev_poly(f,t))
        got[r]=m.add((1,got[r]),(w,P))
assert got==m.coeff

def choose(columns):
    # Row pivots with exact fraction elimination, independent of prime.
    piv={}
    keep=[]
    for src in columns:
        P=src.copy()
        while P:
            w=min(P)
            c=P[w]
            if w not in piv:
                keep.append(w)
                piv[w]={z:v/c for z,v in P.items()}
                break
            for z,a in piv[w].items():
                v=P.get(z,0)-a*c
                if v:P[z]=v
                else:P.pop(z,None)
    return keep

def det(A,p):
    M=[[m.mod(v,p) for v in row] for row in A];out=1
    for i in range(len(M)):
        j=next((j for j in range(i,len(M)) if M[j][i]),None)
        if j is None:return 0
        if j!=i:M[i],M[j]=M[j],M[i];out=-out
        s=M[i][i];out=out*s%p;q=pow(s,-1,p)
        tail=[(k,a) for k,a in enumerate(M[i]) if k>i and a]
        for row in M[i+1:]:
            b=row[i]*q%p
            if b:
                for k,a in tail:row[k]=(row[k]-b*a)%p
                row[i]=0
    return out
cols=[A|B for A,B in zip(m.constraints,m.M)]
constraint_rows=choose(m.constraints)
rows=choose(cols)
assert len(constraint_rows)==290 and len(rows)==len(cols)==300
MM=[[X.get(w,F()) for X in cols] for w in rows]
pp=(1000003,2147483647)
res=[det(MM,p) for p in pp]
assert all(res)
receipt={'bare_prefactors_match_independent_dmk_source':True,'coefficient_terms_by_flavor':list(map(len,got)), 'constraint_rank_Q':len(constraint_rows), 'raw_columns':len(cols),'stacked_rank_Q':len(rows),'primes':pp,'original_minor_determinants':res,'original_rows':rows}
Path(__file__).with_name('cartan-300-minor.json').write_text(json.dumps(receipt,separators=(',',':'))+'\n')
print({k:v for k,v in receipt.items() if k!='original_rows'})
