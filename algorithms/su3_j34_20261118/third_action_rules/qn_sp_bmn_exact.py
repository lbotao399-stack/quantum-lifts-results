"""Path-form Q1,Q2,Q3 on derivative-free inputs, in the pinned code convention.

Source: n4_supercharges.tex eqs. direct, Q1direct, Q2direct, Q3direct.
The first structure constant in eq. direct is reversed to [T_R,C_D].
Pairing a path with its reversal and expressing the physical coefficients
in the existing Q1=1 convention gives (-1)^(n+1) 2^n/((n!)^3(n+1)):
1, -1/6, 1/108. No physical coupling normalization is inferred here.
Only input jets of order zero are implemented.
"""
import sp_environment
from q2d_sp_exact import contract,_reorder_sign
from fortuity.full.words import canon_word,canon_multi
from fortuity.full.letters import is_ferm
from functools import lru_cache
from itertools import product,combinations,permutations
from fractions import Fraction
from math import factorial,comb
from collections import defaultdict


@lru_cache(None)
def theta_rule(inputs):
    n=len(inputs)-1
    if sum(map(len,inputs))<3*n:return ()
    fields=[('F',i) for i in range(n+1)]
    field_par={x:(len(inputs[i])+1)%2 for i,x in enumerate(fields)}
    parity=lambda x:field_par[x] if x in field_par else 1
    current=[x for i,S in enumerate(inputs) for x in ([(i,j) for j in S]+[fields[i]])]
    target=[(i,j) for i,S in enumerate(inputs) for j in S]+fields
    lhs=_reorder_sign(current,target,parity)*(-1)**sum(map(len,inputs))
    out=defaultdict(int)
    for choices in product((0,1),repeat=3*n):
        pref=[(i+choices[3*i+j],j) for i in range(n) for j in range(3)]
        if len(set(pref))!=len(pref):continue
        selected=[{j for i0,j in pref if i0==i} for i in range(n+1)]
        if any(selected[i]!=set(inputs[i]) for i in range(1,n)):continue
        if not selected[0]<=set(inputs[0]) or not selected[n]<=set(inputs[n]):continue
        SD=tuple(j for j in inputs[0] if j not in selected[0])
        SE=tuple(j for j in inputs[n] if j not in selected[n])
        rpar=lambda x:(len(SD)+1)%2 if x=='D' else (len(SE)+1)%2 if x=='E' else 1
        current=pref+[(0,j) for j in SD]+['D']+[(n,j) for j in SE]+['E']
        target=[(i,j) for i,S in enumerate(inputs) for j in S]+['D','E']
        out[SD,SE]+=lhs*(-1)**sum(choices)*_reorder_sign(current,target,rpar)
    return tuple((c,SD,SE) for (SD,SE),c in out.items() if c)


def Qn_bmn_multi(multi,N,n,group='SU',both_paths=False):
    assert n in (1,2,3)
    if group not in ('SU','SO_odd','Sp'):raise NotImplementedError(group)
    letters=[L for w in multi for L in w]
    if any(L[0] or L[1] for L in letters):raise NotImplementedError('Input jets are not implemented')
    count=len(letters);par=[int(is_ferm(L)) for L in letters]
    ids=[];offset=0
    for w in multi:ids.append(list(range(offset,offset+len(w))));offset+=len(w)
    coefficient=Fraction((-1)**(n+1)*2**n,factorial(n)**3*(n+1))
    if both_paths:coefficient/=2
    out=defaultdict(Fraction)
    for chosen in combinations(range(count),n+1):
        rest=[i for i in range(count) if i not in chosen]
        for path in permutations(chosen):
            if not both_paths and path>path[::-1]:continue
            rule=theta_rule(tuple(letters[i][2] for i in path))
            if not rule:continue
            left,right=path[0],path[-1];front=_reorder_sign(list(range(count)),list(path)+rest,lambda i:par[i])
            tokens=['X%d'%i for i in range(n)];reference=[left,right]+rest
            for c0,SD,SE in rule:
                newpar={left:(len(SD)+1)%2,right:(len(SE)+1)%2}
                parity=lambda i:newpar.get(i,par[i])
                for orientation in product((0,1),repeat=n+1):
                    replacements={left:[tokens[0],left],right:[tokens[-1],right]}
                    for i in range(1,n):replacements[path[i]]=[tokens[i-1],tokens[i]]
                    for i,t in enumerate(path):
                        if orientation[i]:replacements[t]=replacements[t][::-1]
                    necks=[[v for t in word for v in replacements.get(t,[t])] for word in ids]
                    terms=[(1,necks)]
                    for token in tokens:
                        terms=[(s*sg,new) for s,old in terms for sg,new in contract(old,token,group)]
                    for sign,necks in terms:
                        if any(len(t)==1 for t in necks):continue
                        factor=N**sum(not t for t in necks);necks=[t for t in necks if t]
                        written=[t for w in necks for t in w]
                        base=coefficient*c0*front*(-1)**sum(orientation)*sign*factor*_reorder_sign(written,reference,parity)
                        for k in range(n+1):
                            LD=(n-k,k,SD);LE=(k,n-k,SE)
                            cf=base*(-1)**k*comb(n,k)*(factorial(k)*factorial(n-k))**2;ws=[]
                            for word in necks:
                                sg,w=canon_word(tuple(LD if i==left else LE if i==right else letters[i] for i in word))
                                if not sg:cf=0;break
                                cf*=sg;ws.append(w)
                            if cf:
                                sg,w=canon_multi(ws)
                                if sg:out[w]+=cf*sg
    return {w:c for w,c in out.items() if c}
