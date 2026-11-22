"""Exact SU3→SU2 block classical separating witness for public R007."""
from fractions import Fraction as F
from collections import defaultdict
from functools import lru_cache
from itertools import product,combinations
from hashlib import sha256
from pathlib import Path
import json,gzip
import sys
assert len(sys.argv)==3, 'provide frozen O32_C_even.json.gz and o0_classical_certificate.json'
src,dual=map(Path,sys.argv[1:])
assert sha256(src.read_bytes()).hexdigest()=='40c4cbb876c01811d0dd83035ab9080ad71ebb603331304b309f4d0d0335e623'
assert sha256(dual.read_bytes()).hexdigest()=='78c4eab761234f9c136a42c1481a3c92cc03887e8345545c0e6810995e019eda'
def canon(xs):
    odd=[x for x in xs if x>=12]
    if len(set(odd))<len(odd):return (),0
    s=(-1)**sum(odd[j]>odd[k] for k in range(len(odd)) for j in range(k))
    return tuple(sorted(xs)),s
@lru_cache(None)
def tr(w):
    groups=[defaultdict(F) for i in range(5)]
    zeros=[i for i,x in enumerate(w) if x==0]
    for k in range(min(4,len(zeros))+1):
        for removed in combinations(zeros,k):
            keep=[x for i,x in enumerate(w) if i not in removed]
            if not keep:
                groups[k][()]+=2+(-2)**k
                continue
            for colors in product(range(3), repeat=len(keep)):
                last=-1; phase=0
                for c in colors:
                    if last<0:last=c;continue
                    if last==c:last=-1;continue
                    phase+=1 if c==(last+1)%3 else -1
                    last=3-last-c
                if last>=0:continue
                key,sgn=canon([3*x+c if x<3 else 12+3*(x-3)+c for x,c in zip(keep,colors)])
                if sgn:
                    # phase of Pauli trace: 2 * i**phase. Track the power
                    # separately until pairing; do not use complex floats.
                    groups[k][(key,phase%4)]+=2*sgn
    for k in range(5):
        # For pure shift traces use the same data form as all others.
        if () in groups[k]:groups[k][((),0)]+=groups[k].pop(())
    return groups
def mul(A,B):
    P=[defaultdict(F) for i in range(5)]
    for j in range(5):
        for k in range(5-j):
            for (a,p),c in A[j].items():
                for (b,q),d in B[k].items():
                    z,s=canon(a+b)
                    if s:P[j+k][(z,(p+q)%4)]+=c*d*s
    return [{w:c for w,c in col.items() if c} for col in P]
def w_of(x):
    e,m=x
    return tuple([i for i,k in enumerate(e) for _ in range(k)]+[12+i for i in range(9) if m>>i&1])
C=json.loads(dual.read_bytes())
ell={w_of(C['target_basis'][i]):F(c) for i,c in C['dual']}
O0={w_of(C['target_basis'][i]):F(c) for i,c in C['O0']}
P=defaultdict(F)
for c,m in json.loads(gzip.decompress(src.read_bytes()))['terms']:
    ms=[{((),0):F(c)},{},{},{},{}]
    for w in m:ms=mul(ms,tr(tuple(w)))
    for (w,k),a in ms[4].items():P[w,k]+=a
part=[defaultdict(F),defaultdict(F)]
for (w,k),c in P.items():part[k%2][w]+=c*(-1)**(k//2)
part=[{w:c for w,c in p.items() if c} for p in part]
pairs=[sum(c*ell.get(w,0) for w,c in p.items()) for p in part]
assert any(pairs)
match=[]
for p in part:
    if not p:match.append('0');continue
    ratios={p.get(k,0)/v for k,v in O0.items()}
    match.append(str(next(iter(ratios))) if p.keys()==O0.keys() and len(ratios)==1 else None)
assert pairs==[0,244944000] and match==['0','-8164800']
print(json.dumps({'source_sha256':sha256(src.read_bytes()).hexdigest(),'classical_dual_sha256':sha256(dual.read_bytes()).hexdigest(),'shift':'diag(X + h I_2, -2h) on phi1 only, coefficients of h^4','generator':'unscaled Pauli matrices','real_imag_support':[len(x) for x in part],'real_imag_pairing':[str(x) for x in pairs],'literal_multiple_of_saved_O0':match,'literal_O0_pairing':str(sum(v*ell.get(k,0) for k,v in O0.items()))},indent=2))
