"""The scalar--matter singlet and the seven-fermion classical correction to its fifth-power action."""
from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path
import argparse

import numpy as np

import mixed_incoming as m
from f4_colour import basis
from rank_four_s5 import decode_rows
from q1d_exact import pair_terms

G=((1,0,()),(0,1,()))


def singlet():
    result=Counter()
    for i in range(3):
        scalar=(0,0,(i,))
        matter=(0,0,tuple(j for j in range(3) if i!=j))
        s,w=m.so_canon(((scalar,matter),))
        result[w]+=s*(-1)**i
    for i in (1,2,3):
        action=Counter()
        for w,c in result.items():
            for v,n in m.raise_word(w,i).items():
                action[v]+=c*n
        assert all(n==0 for n in action.values())
    return dict(result)


def multiply(rows, power):
    out=Counter()
    for w,c in rows.items():
        sign,word=m.so_canon(w+(G,)*power)
        if sign:
            out[word]+=c*sign
    return {w:c for w,c in out.items() if c}


def verify():
    generators,_=basis()
    gram=np.einsum("iab,jba->ij",generators,generators)
    assert np.array_equal(gram,np.diag(gram.diagonal()))
    den=lcm(*(abs(int(n)) for n in gram.diagonal()))
    dual=generators*(den//gram.diagonal())[:,None,None]
    assert np.array_equal(np.einsum("iab,jba->ij",generators,dual),den*np.eye(52,dtype=int))
    for x in generators:
        comm=dual@x-x@dual
        assert np.array_equal((generators@comm-comm@generators).sum(axis=0),3*den*x)
    source=singlet()
    # Two exterior generators isolate the coefficient of the ordered gaugino pair.
    # An exchanged odd output reverses that coefficient. The trace of two
    # commutators divided by the dual denominator is computed directly.
    x=generators[0]
    one=generators@x-x@generators
    two=dual@x-x@dual
    result=Fraction(0)
    for (word,),c in source.items():
        for coeff,left,right in pair_terms(*word):
            assert (left,right) in (G,G[::-1])
            exterior=1 if (left,right)==G else -1
            result+=c*coeff*exterior*Fraction(int(np.einsum("iab,iba->",one,two)),den*int(np.trace(x@x)))
    assert result==18
    weighted=multiply(source,4)
    assert all(not pair_terms(l,g) and not pair_terms(g,l) for w in weighted for tr in w for l in tr for g in G)
    variation=m.q(weighted)
    assert variation==multiply(m.q(source),4) and len(variation)==1
    return dict(source=[[m.encode(w),c] for w,c in sorted(weighted.items())],
                classical_image=[[m.encode(w),c] for w,c in sorted(variation.items())],
                colour_dual_denominator=den,adjoint_casimir=3,casimir_checks=52,
                first_action_fifth_power=int(result),field="Q")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    assert __debug__
    m.write(args.out,verify())
