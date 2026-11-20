"""Pinned Q3 path action extended to input derivatives by exact kernel moments.

Colour reconnections, theta rules, Koszul signs and the coefficient 1/108
are copied unchanged from qn_sp_bmn_exact.Qn_bmn_multi at n=3.
"""
import sp_environment
from q2d_sp_exact import contract,_reorder_sign
from qn_sp_bmn_exact import theta_rule
from q3_jet_moments import jet_terms
from fortuity.full.words import canon_word,canon_multi
from fortuity.full.letters import is_ferm
from functools import lru_cache
from itertools import product,combinations,permutations
from fractions import Fraction
from math import factorial,comb
from collections import defaultdict


def Q3d_multi(multi,N,group='Sp',both_paths=False):
    n=3
    if group not in ('SU','SO_odd','Sp'):raise NotImplementedError(group)
    letters=[L for w in multi for L in w]
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
            jets=jet_terms(tuple(letters[i][:2] for i in path))
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
                        for jet_coefficient, first_jets, second_jets in jets:
                            LD=first_jets+(SD,);LE=second_jets+(SE,)
                            cf=base*jet_coefficient;ws=[]
                            for word in necks:
                                sg,w=canon_word(tuple(LD if i==left else LE if i==right else letters[i] for i in word))
                                if not sg:cf=0;break
                                cf*=sg;ws.append(w)
                            if cf:
                                sg,w=canon_multi(ws)
                                if sg:out[w]+=cf*sg
    return {w:c for w,c in out.items() if c}
