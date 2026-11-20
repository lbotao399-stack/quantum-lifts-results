"""Raise the known eight-letter SO(8) BMN class and return it to derivative grade zero."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import chain
import json
from pathlib import Path

import mixed_incoming as m
from fortuity.full.letters import Q_letter, is_ferm
from gaugino_primary_order import STEPS, charge


@lru_cache(None)
def substitution(letter, spin, flavour, increase):
    i,j,theta = letter
    sign = (-1)**sum(x > flavour for x in theta)
    if increase:
        if flavour not in theta:
            return 0, None
        return sign*((i,j)[spin]+1), (i+(spin==0),j+(spin==1),tuple(x for x in theta if x!=flavour))
    if not (i,j)[spin] or flavour in theta:
        return 0, None
    return sign, (i-(spin==0),j-(spin==1),tuple(sorted((*theta,flavour))))


def check_letter(letter,spin,flavour,increase):
    n,changed=substitution(letter,spin,flavour,increase)
    terms=Counter()
    if n:
        for c,left,right in Q_letter(changed):
            terms[left,right]+=c*n
    for c,left,right in Q_letter(letter):
        a,new=substitution(left,spin,flavour,increase)
        if a: terms[new,right]+=a*c
        b,new=substitution(right,spin,flavour,increase)
        if b: terms[left,new]+=(-1)**is_ferm(left)*b*c
    assert not any(terms.values()), (letter,spin,flavour,increase)


def substitute(expr, spin, flavour, increase):
    out=Counter()
    for word,value in expr.items():
        phase=1
        for i,tr in enumerate(word):
            for j,letter in enumerate(tr):
                n,new=substitution(letter,spin,flavour,increase)
                if n:
                    sg,converted=m.so_canon(word[:i]+(tr[:j]+(new,)+tr[j+1:],)+word[i+1:])
                    if sg: out[converted]+=value*phase*sg*n
                if is_ferm(letter): phase=-phase
    return {w:v for w,v in out.items() if v}


def run(path):
    data=json.loads(gzip.decompress(path.read_bytes()))
    block=data["bases"]["n1_r0"]
    words=[m.decode(w) for w in block["words"]]
    co=data["representatives"]["O18_plus_g2_O1"]["coefficients"]
    original=Counter()
    for scale,row in zip(co[:len(block["rows"])],block["rows"],strict=True):
        for j,n in row: original[words[j]]+=n*Fraction(scale)
    assert all(v.denominator==1 for v in original.values())
    original={w:int(v) for w,v in original.items() if v}
    assert len(original)==579
    expr=original
    records=[]
    for increase,step in chain(((True,s) for s in STEPS),((False,s) for s in reversed(STEPS))):
        alphabet={l for w in expr for tr in w for l in tr}
        for l in alphabet: check_letter(l,*step,increase)
        expr=substitute(expr,*step,increase)
        if not expr:
            records.append(dict(increase=increase,step=step,terms=0,letter_checks=len(alphabet)))
            return dict(input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),steps=records,return_coefficient=0)
        weight,length=charge(expr)
        records.append(dict(increase=increase,step=step,terms=len(expr),weight=weight,length=length,
                            letter_checks=len(alphabet)))
        print(increase,step,len(expr),flush=True)
        if len(records)==6:
            assert weight==(3,3,1,1,1) and length==8
            assert all(not any(sum_images.values()) for sum_images in raising_images(expr))
            digest=hashlib.sha256(m.compact([[m.encode(w),n] for w,n in sorted(expr.items())])).hexdigest()
    first=min(original)
    coefficient=Fraction(expr[first],original[first])
    assert coefficient and set(expr)==set(original)
    assert all(n*coefficient==expr[w] for w,n in original.items())
    return dict(input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),steps=records,
                target_sha256=digest,return_coefficient=str(coefficient))


def raising_images(expr):
    for axis in (1,2,3):
        out=Counter()
        for w,n in expr.items():
            for term,c in m.raise_word(w,axis).items(): out[term]+=n*c
        yield out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",type=Path,default=Path(__file__).with_name("so8_level18")/"perturbative_complex.json.gz")
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    assert __debug__
    record=run(args.input)
    m.write(args.out,record)
    print("return coefficient",record["return_coefficient"])
