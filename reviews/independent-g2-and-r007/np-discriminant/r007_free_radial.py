"""Free-singleton lowering and Cartan Gram tests on Dot's exact SU(3) outgoing state."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from math import factorial, prod
from pathlib import Path
import gzip
import json

import r007_primary_replay as jets


def source_input(name, digest):
    data=Path(__file__).with_name(name).read_bytes()
    assert sha256(data).hexdigest()==digest
    return data


def target():
    recipes=json.loads(gzip.decompress(source_input('su30_0.json.gz',
        '81f936300745c652c4279bb769365b59df01fdcb321860d657b7b6d3e1188ebc')))['selected_recipes']
    selected=json.loads(source_input('R007_recipe.json',
        'c82ec0745328b8482c12dec843bb25d25345f9f7dffa6776720a37971acba4df'))['source_linear_combination']
    all_terms=Counter()
    for i,c in selected:
        term={():Fraction(c)}
        for g in recipes[i]:
            term=jets.multiply(term,jets.factor(tuple(g)))
        all_terms.update(term)
    out=jets.clean(all_terms)
    sizes=[len(out)]
    for dot,flavour in ((0,1),(1,1),(0,2),(1,2)):
        out=jets.act(out,('T',dot,flavour))
        sizes.append(len(out))
    assert sizes==[348,1059,1864,4561,5080]
    return out


def oscillator(letter):
    p,q,holes=letter
    num_a=p+q+len(holes)-1
    assert num_a>=0
    occupied=tuple(i for i in range(3) if i not in holes)
    assert num_a-p-q+len(occupied)==2
    sign=(-1)**(sum(holes)+len(holes)*(len(holes)-1)//2)
    return (num_a,p,q,occupied),Fraction(sign,factorial(p)*factorial(q))


def generator(letter, op):
    state,coefficient=oscillator(letter)
    num_a,p,q,occupied=state
    occupation=[num_a,p,q]
    filled=set(occupied)
    kind,i,j=op
    if kind=='E': actions=[('c',i),('c+',j)]; coefficient*=-1
    elif kind=='J': actions=[('b',1),('b+',0)]
    elif kind=='K': actions=[('b',i),('a',0)]
    elif kind=='S': actions=[('a',0),('c+',i)]; coefficient*=-1
    elif kind=='T': actions=[('b',i),('c',j)]
    elif kind=='D' and i<2: actions=[('b+',i),('a+',0)]
    elif kind=='D': actions=[('c',i-2),('a+',0)]
    else: raise ValueError(op)
    for symbol,k in actions:
        if symbol[0]=='c':
            if (k in filled)!=(symbol=='c'):
                return ()
            coefficient*=(-1)**sum(h<k for h in filled)
            if symbol=='c': filled.remove(k)
            else: filled.add(k)
        else:
            axis=0 if symbol[0]=='a' else k+1
            if symbol.endswith('+'):
                occupation[axis]+=1
            else:
                coefficient*=occupation[axis]
                if not coefficient: return ()
                occupation[axis]-=1
    next_letter=occupation[1],occupation[2],tuple(i for i in range(3) if i not in filled)
    check,change=oscillator(next_letter)
    assert check==(occupation[0],occupation[1],occupation[2],tuple(sorted(filled)))
    return ((coefficient/change,(next_letter,)),)


def inner_letter(letter):
    (num_a,p,q,_),coefficient=oscillator(letter)
    return coefficient**2*factorial(num_a)*factorial(p)*factorial(q)


@lru_cache(None)
def orthogonal_trace(trace):
    # H1=diag(1,-1,0), H2=diag(1,1,-2). The defining trace form is diag(2,6).
    symbols=tuple((p,q,sum(1<<i for i in holes)) for p,q,holes in trace)
    out=Counter()
    for colours in product(range(2),repeat=len(trace)):
        n=colours.count(0)
        tr=1+(-1)**n+(int(n==0)*(-2)**len(trace))
        if not tr: continue
        a=tuple(x for x,k in zip(symbols,colours) if k==0)
        b=tuple(x for x,k in zip(symbols,colours) if k==1)
        sign1,first=jets.ordered(a,jets.diagonality)
        sign2,second=jets.ordered(b,jets.diagonality)
        crossing=sum(jets.diagonality(symbols[i])*jets.diagonality(symbols[j])
                     for i,j in combinations(range(len(trace)),2) if colours[i]>colours[j])
        out[first,second]+=sign1*sign2*(-1)**crossing*tr
    return jets.clean(out)


def inner_polynomial(poly):
    total=Fraction()
    for (first,second),coeff in poly.items():
        term=coeff**2
        for group,colour_length in ((first,2),(second,6)):
            for (p,q,mask),n in Counter(group).items():
                holes=tuple(i for i in range(3) if mask>>i&1)
                assert n==1 or not jets.odd((p,q,holes))
                term*=factorial(n)*(inner_letter((p,q,holes))/colour_length)**n
        total+=term
    return total


def project_orthogonal(state):
    old=jets.diag
    try:
        jets.diag=orthogonal_trace
        return jets.rho(state)
    finally:
        jets.diag=old


def run():
    ops=[('E',0,1),('E',1,2),('J',0,0)]+[('T',d,i) for d,i in product(range(2),range(3))]
    ops += [('S',i,0) for i in range(3)]+[('K',d,0) for d in range(2)]
    translations=[('D',i,0) for i in range(5)]
    # Includes the entire alphabet of the pinned polynomial and its elementary insertions.
    letters=[(p,q,i) for p,q in product(range(5),repeat=2) for k in range(4)
             for i in combinations(range(3),k) if 0<p+q+k<=7]
    for letter in letters:
        for op in ops+translations:
            assert generator(letter,op)==jets.operation(letter,op),(letter,op)
    source=target()
    assert {l for traces in source for tr in traces for l in tr}<=set(letters)
    assert all(not jets.act(source,op) for op in ops)
    assert not jets.act(source,('Q',0,0))
    ortho=project_orthogonal(source)
    gram=inner_polynomial(ortho)
    assert gram>0
    data=[[w,str(c)] for w,c in sorted(ortho.items())]
    digest=sha256(json.dumps(data,separators=(',',':')).encode()).hexdigest()
    return dict(source_terms=len(source),tested_letters=len(letters),operations=len(ops+translations),
                singleton_equations=len(letters)*len(ops+translations),target_lowering_equations=len(ops),
                Q0_closed=True,orthogonal_Cartan_terms=len(ortho),orthogonal_Cartan_sha256=digest,
                Cartan_Gram_defining_trace_unit_singleton=str(gram))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args()
    if not __debug__: raise RuntimeError('run with assertions enabled')
    receipt=run()
    args.out.write_text(json.dumps(receipt,indent=2)+'\n')
    print(receipt,flush=True)
