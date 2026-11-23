#!/usr/bin/env python3
"""Closed two-coefficient dual in the fixed native scalar Pfaffian convention.

The proof and its scope are in UNIVERSAL_TWO_PROBE_DUAL.md. --verify compares
the formula against actual ordered-trace coefficient extraction and actual Q0.
"""
import argparse
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
import json
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'scalar'))
import minimal_scalar_certificates as minimal
import pfaffian_scalar_candidates as p


def double_factorial(n):
    assert n >= -1
    return prod(range(n,0,-2))


def terminal_coefficient(shape, triples):
    x,y,z=shape
    assert x>=y>=z>=3 and all(v%2 for v in shape)
    assert 1<=triples<=z and triples%2
    value=Fraction((-1)**((x-z)//2)*double_factorial(x+2)*double_factorial(y+1)*double_factorial(z),4**((z-1)//2))
    for r in range(1,triples,2):
        value*=Fraction(4*(x-r+2)*(y-r+1)*(z-r),(r+1)*(r+2))
    return value


def word(shape, triples):
    raw=[(0,1,2)]*triples
    raw += [(i,i) for i,n in enumerate(shape) for _ in range((n-triples)//2)]
    sign,canonical=p.canonical(raw)
    assert sign==1
    return canonical


def run(shape,verify):
    x,y,z=shape
    values=[terminal_coefficient(shape,r) for r in (1,3)]
    dual=[Fraction(3,2)/values[0],-Fraction(1,2)/values[1]]
    divergences=[sum(shape)+6]+([x+4] if x>y else [])+([z] if y>z else [])
    denominator=z*(y+1)*(x+2)
    matrix=[[Fraction(-2*r*d,denominator)*value for d in divergences]+[value] for r,value in zip((1,3),values)]
    assert [sum(dual[i]*matrix[i][j] for i in range(2)) for j in range(len(divergences)+1)]==[0]*len(divergences)+[1]
    checked=[]
    if verify:
        sources,target,_,_=minimal.operators(shape)
        sources=list(map(p.Coefficient,sources));target=p.Coefficient(target)
        for r in range(1,z+1,2):
            w=word(shape,r)
            value=terminal_coefficient(shape,r)
            assert target(w)==value,(shape,r,'target')
            actual=[p.image(source,w) for source in sources]
            expected=[Fraction(-2*r*d,denominator)*value for d in divergences]
            assert actual==expected,(shape,r,'image')
            checked.append(r)
    return {'shape':shape,'SO_colours':sum(shape)-2,'length':sum(shape)-1,'j':2*sum(shape),
            'ordered_probe_words':[p.encode(word(shape,r)) for r in (1,3)],
            'two_probe_matrix':[[str(v) for v in row] for row in matrix],
            'dual_coefficients':list(map(str,dual)),
            'pairings':[0]*len(divergences)+[1],
            'actual_coefficient_checks_at_ABC_trace_counts':checked}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--shape',type=int,nargs=3,action='append');ap.add_argument('--verify',action='store_true');ap.add_argument('--out',type=Path,default=HERE/'CLOSED_DUAL_VERIFIED.json');a=ap.parse_args()
    started=time.monotonic()
    shapes=a.shape or [(3,3,3),(5,3,3),(7,3,3),(9,3,3),(11,3,3),(5,5,3),(7,5,3),(9,5,3),(7,7,3),(5,5,5),(7,5,5),(9,5,5),(7,7,5),(7,7,7),(9,9,3),(11,7,5),(9,7,7),(9,9,9),(13,5,3)]
    records=[run(tuple(shape),a.verify) for shape in shapes]
    result={'passed':True,'arithmetic':'Q','physical_normalization':'unchanged native Q0 and auxiliary H,E,F conventions','records':records,'seconds':time.monotonic()-started}
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'shapes':len(records),'actual_trace_counts_checked':sum(len(r['actual_coefficient_checks_at_ABC_trace_counts']) for r in records),'seconds':result['seconds']},indent=2))


if __name__=='__main__': main()
