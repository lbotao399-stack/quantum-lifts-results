"""Four Cartan probes on the ten-letter singlets and on products of single gravitons.

Use Chang--Yin's graded-symmetric single traces of first derivatives of Psi,
with at most one outside left theta derivative. All flavour and spin indices
are contracted with the original alternating tensors. The colour checks use
the defining traces in SO(7) and USp(6).
"""
import argparse
import json
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path

from sympy import Matrix
import mixed_incoming as m
from mixed33_exact import key, read
from fortuity.full.letters import Q_letter, is_ferm, letters_up_to

AXES = ((1,0), (0,1))
# The three entries count inside theta derivatives, inside z derivatives,
# and outside theta derivatives, respectively. Each inside derivative acts
# on its own superfield. An outside derivative acts on the entire trace.
U, L, G, V = (2,0,0), (1,1,0), (1,1,1), (2,0,1)
FAMILIES = ((U,U,U,L,G), (U,U,V,L,L), (U,U,V,(2,2,0)), (U,(3,1,0),(3,1,1)))
PROBES = ((0,0,1,1), (0,1,0,1), (0,0,0,0), (None,None,0,None))


def sign(seq):
    return (-1)**sum(i>j for at,i in enumerate(seq) for j in seq[at+1:])


def wedge(left,right):
    a,i = left
    b,j = right
    if i & j:
        return 0,0
    swaps = sum((j & ((1 << k)-1)).bit_count() for k in range(i.bit_length()) if (i>>k)&1)
    return a*b*(-1)**swaps, i|j


def fields(probe):
    r,u,i,j = probe
    jet = {(0,0,(f,)): (f,1,0) for f in range(3)}
    jet[1,0,()] = i,1,1
    if r is None:
        jet[0,1,(2,)] = 1,1,0
        return jet,1
    jet[0,1,()] = j,1,2
    jet[0,0,tuple(f for f in range(3) if f != u)] = r,(-1)**u,4
    return jet,7


def trace_value(tr, jet, group):
    if len(tr)%2 or any(l not in jet for l in tr):
        return 0,0
    # Each H_i is supported on its own two-dimensional block. There it is
    # [[0,1],[-1,0]] in SO(7) or diag(1,-1) in USp(6).
    if len({jet[l][0] for l in tr}) != 1:
        return 0,0
    out = (2*(-1)**(len(tr)//2) if group == "SO(7)" else 2),0
    for l in tr:
        out = wedge(out,jet[l][1:])
    return out


def expression_value(expr, jet, group):
    row = Counter()
    for w,c in expr.items():
        term = c,0
        for tr in w:
            term = wedge(term,trace_value(tr,jet,group))
        row[term[1]] += term[0]
    return {i:c for i,c in row.items() if c}


def derivative(word, flavour):
    if flavour is None:
        yield 1,word
        return
    phase = 1
    for i,(a,b,theta) in enumerate(word):
        if flavour not in theta:
            # d_flavour d_theta Psi on a scalar first derivative.
            factor = (-1 if theta and theta[0] > flavour else 1)*phase
            yield factor, word[:i]+((a,b,tuple(sorted(theta+(flavour,)))),)+word[i+1:]
        if is_ferm(word[i]):
            phase = -phase


@lru_cache(None)
def factor(spec, flavours, spins):
    n,h,d = spec
    inside = tuple((0,0,(i,)) for i in flavours[:n]) + tuple((*AXES[j],()) for j in spins)
    outside = flavours[-1] if d else None
    word_sum = Counter()
    for perm in permutations(range(n+h)):
        phase = sign([i for i in perm if i >= n])
        for sg,tr in derivative(tuple(inside[i] for i in perm),outside):
            cyclic,w = m.so_canon((tr,))
            if cyclic:
                word_sum[w] += phase*sg*cyclic
    expr = {w:Fraction(c,math.factorial(n+h)) for w,c in word_sum.items() if c}
    assert not m.q(expr)
    return expr


def components(spec, spins, jet, group):
    n,h,d = spec
    assert h == len(spins)
    out = []
    for fl in product(range(3),repeat=n+d):
        expr = factor(spec,fl,spins)
        v = expression_value(expr,jet,group)
        # On these commuting Cartan elements the graded symmetrization can
        # also be evaluated in its original order, before differentiation.
        inside = tuple((0,0,(i,)) for i in fl[:n])+tuple((*AXES[j],()) for j in spins)
        control = Counter()
        for sg,tr in derivative(inside,fl[-1] if d else None):
            c,mask = trace_value(tr,jet,group)
            control[mask] += sg*c
        assert v == {a:b for a,b in control.items() if b}
        out.extend((fl,(c,mask)) for mask,c in v.items())
    return out


def tensor(specs, probe, group):
    jet,mask = fields(probe)
    values = Counter()
    for spins,sg in (((0,1),1),((1,0),-1)):
        p = 0
        rows = []
        for spec in specs:
            h = spec[1]
            rows.append(components(spec,spins[p:p+h],jet,group))
            p += h
        assert p == 2
        for pieces in product(*rows):
            term = sg,0
            labels = ()
            for f,v in pieces:
                term = wedge(term,v)
                labels += f
            if term[1] == mask and term[0]:
                assert len(labels) == 9
                values[labels] += term[0]
    return {k:v for k,v in values.items() if v}


def partitions(ids):
    if not ids:
        yield ()
        return
    for i,j in combinations(ids[1:],2):
        tr = (ids[0],i,j)
        for other in partitions(tuple(k for k in ids if k not in tr)):
            yield (tr,)+other


def contract(row, triples):
    result = Fraction(0)
    for f,v in row.items():
        for tri in triples:
            labels = [f[i] for i in tri]
            if len(set(labels)) != 3:
                break
            v *= sign(labels)
        else:
            result += v
    return result


def cartan_differentials():
    # Tensoring commuting Cartan generators with exterior coefficients makes
    # every graded commutator vanish. Check the ordered Q0 letter rules too.
    checked = 0
    for probe in PROBES:
        jet,_ = fields(probe)
        for l in letters_up_to(1,1):
            values = Counter()
            for c,a,b in Q_letter(l):
                if a in jet and b in jet and jet[a][0] == jet[b][0]:
                    v,mask = wedge(jet[a][1:],jet[b][1:])
                    values[jet[a][0],mask] += c*v
            assert not any(values.values())
            checked += 1
    return checked


def run(formal):
    z = read(formal)
    assert z["mu"] == [1,1,3,3,3] and z["length"] == 10
    cycles = []
    for b in z["blocks"]:
        for rows in b["representatives"]:
            expr = {key(w):c for w,c in rows}
            assert not m.q(expr)
            cycles.append(expr)
    tests = tuple(partitions(tuple(range(9))))
    assert len(tests) == 280
    results = {}
    for group in ("SO(7)","USp(6)"):
        classical = [[expression_value(expr,jet,group).get(mask,0) for jet,mask in map(fields,PROBES)] for expr in cycles]
        assert Matrix(classical).rank() == 4
        matrix, contractions = [], []
        for family in FAMILIES:
            values = [tensor(family,probe,group) for probe in PROBES]
            for tri in tests:
                row = [contract(v,tri) for v in values]
                if Matrix(matrix+[row]).rank() > len(matrix):
                    matrix.append(row)
                    contractions.append(tri)
                    break
            else:
                raise ValueError(("family failed to add a row",family,matrix))
        results[group] = dict(formal_class_evaluations=classical,
                              formal_minor_indices=[0,1,2,6],
                              formal_minor=str(Matrix([classical[i] for i in (0,1,2,6)]).det()),
                              product_matrix=[list(map(str,row)) for row in matrix],
                              product_determinant=str(Matrix(matrix).det()),epsilon_triples=contractions)
    return dict(weight=z["mu"],length=z["length"],four_cartan_probes=PROBES,
                Q0_letter_checks=cartan_differentials(),signed_trace_terms=list(map(len,cycles)),
                factors_by_inside_theta_z_and_outside_theta=FAMILIES,groups=results)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--formal",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    assert __debug__
    cert=run(a.formal)
    a.out.write_text(json.dumps(cert,indent=2)+"\n")
    print({g:r["product_determinant"] for g,r in cert["groups"].items()})
