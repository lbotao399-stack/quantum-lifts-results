"""Test a complete scalar layer on arbitrary three-colour matrices commuting with the selected scalar."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import lcm
from pathlib import Path

from su3_level32_layers import DATA, content, formal, joint, residual
from su3_level32_levi import TRANSLATION
from su3_level32_trace_factors import evaluator
from verify import exterior_sign, multiply


@lru_cache(maxsize=30000)
def trace_terms(trace):
    # At diag(1,1,-2), each commuting traceless matrix is diag(A,-tr A).
    active = bytes(k for k in trace if k)
    power = len(trace)-len(active)
    if not active:
        return [((),2+(-2)**power)]
    sign,upper = (1,active) if len(active)==1 else formal.canon(active)
    lower = tuple(bytes((k,)) for k in active)
    reorder,lower = formal.multi(lower)
    rows = Counter()
    if sign:
        rows[(upper,)] += sign
    if reorder:
        rows[lower] += reorder*(-1)**len(active)*(-2)**power
    return [(w,c) for w,c in rows.items() if c]


def restrict_words(words, denominator):
    result = Counter()
    for w,coefficient in words.items():
        integer = coefficient*denominator
        assert integer.denominator==1
        for choice in product(*(trace_terms(t) for t in w)):
            factors = tuple(t for row,_ in choice for t in row)
            sign,canon = formal.multi(factors)
            if sign:
                value = sign*int(integer)
                for _,n in choice:
                    value *= n
                result[canon] += value
    return {w:c for w,c in result.items() if c}


@lru_cache(maxsize=30000)
def trace_polynomial(trace):
    result = Counter()
    for walk in product(range(2),repeat=len(trace)):
        power = mask = 0
        sign = 1
        for k,i,j in zip(trace,walk,walk[1:]+walk[:1],strict=True):
            slot = 4*k+2*i+j
            if formal.ODD[k]:
                bit = 1 << slot
                sign *= exterior_sign(mask,bit)
                if not sign:
                    break
                mask |= bit
            else:
                power += 1 << (4*slot)
        if sign:
            result[power,mask] += sign
    return {k:v for k,v in result.items() if v}


def evaluate(words):
    result = Counter()
    for word,c in words.items():
        row = {(0,0):c}
        for trace in word:
            row = multiply(row,trace_polynomial(trace))
        result.update(row)
    return {k:v for k,v in result.items() if v}


def general_colour_test(target, row, seed):
    diagonal = ((1,0,0,0,0,0,0,0,-1),(0,0,0,0,1,0,0,0,-1))
    basis = diagonal+tuple(tuple(int(k==i) for k in range(9)) for i in (1,2,3,5,6,7))
    ev,assignment = evaluator(target.translate(TRANSLATION),basis,seed)
    value = sum((c*ev(tuple(tr.translate(TRANSLATION) for tr in w)) for w,c in row.items()),Fraction())
    return str(value),assignment


def check_layer(data, correction, power, seed, out):
    saved = joint.s.read(correction)
    terms = {joint.s.unpack_hex(w):Fraction(c) for w,c in saved}
    assert len(terms)==len(saved) and all(content(w).count(0)==5 for w in terms)
    rows = residual(data,terms,lambda t:t.count(0)==power)
    report = {}
    first = None
    for t,row in rows.items():
        den = lcm(*(c.denominator for c in row.values()))
        restricted = restrict_words(row,den)
        polynomial = evaluate(restricted)
        record=dict(three_colour_words=len(row),centralizer_trace_words=len(restricted),
                    nonzero_centralizer_monomials=len(polynomial))
        if polynomial:
            k=min(polynomial)
            record["first_coefficient"]=[str(k[0]),str(k[1]),str(Fraction(polynomial[k],den))]
        if row and first is None:
            value,assignment = general_colour_test(t,row,seed)
            if Fraction(value):
                first=dict(content=t.hex(),coefficient=value,seed=seed,assignment=assignment)
        report[t.hex()]=record
    summary = dict(scalar_power=power,contents=len(rows),formal_zero=sum(not row for row in rows.values()),
                   centralizer_nonzero_contents=sum(bool(r["nonzero_centralizer_monomials"]) for r in report.values()),
                   unrestricted_witness=first,checks=report)
    joint.s.write(out,summary)
    print(power,len(rows),summary["centralizer_nonzero_contents"],first and first["coefficient"],flush=True)


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--correction",required=True,type=Path)
    p.add_argument("--power",type=int,choices=range(3,8),default=5)
    p.add_argument("--seed",type=int,default=894611)
    p.add_argument("--out",required=True,type=Path)
    args = p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    check_layer(args.data,args.correction,args.power,args.seed,args.out)
