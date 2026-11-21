"""Spin and flavour highest weights in the first-differential images of gaugino powers."""
import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import mixed_incoming as m
from fortuity.full.letters import Q_letter, is_ferm, letters_up_to
from gaugino_primary_order import apply, charge, reverse
from rank_four_lowered_pair import check_classical
from rank_four_lowered_replay import derivative
from rank_four_sixth import G


def reverse_letter(letter,spin,flavour):
    x,y,theta=letter
    if flavour not in theta:
        return None
    sign=(-1)**sum(i>flavour for i in theta)
    return sign*(letter[spin]+1),(x+int(spin==0),y+int(spin==1),tuple(i for i in theta if i!=flavour))


def check_reverse(alphabet):
    checks=0
    for letter in alphabet:
        for spin in range(2):
            for flavour in range(3):
                value=Counter()
                lowered=reverse_letter(letter,spin,flavour)
                if lowered:
                    n,new=lowered
                    for c,a,b in Q_letter(new):
                        value[a,b]+=c*n
                for n,a,b in Q_letter(letter):
                    for i,part in enumerate((a,b)):
                        raised=reverse_letter(part,spin,flavour)
                        if raised:
                            c,new=raised
                            key=(new,b) if i==0 else (a,new)
                            value[key]+=n*c*(-1 if i and is_ferm(a) else 1)
                assert not any(value.values()),(letter,spin,flavour,value)
                checks+=1
    return checks


def check_highest(expr):
    for root in (1,2,3):
        value=Counter()
        for word,n in expr.items():
            for target,c in m.raise_word(word,root).items():
                value[target]+=c*n
        assert not any(value.values())


def records(power):
    orientation,w=m.so_canon((G,)*power)
    initial={w:orientation}
    result=[]
    for ascending in combinations_with_replacement(range(3),3):
        mu=ascending[::-1]
        # Replacing minus before plus leaves the highest dotted-spin weight.
        steps=[(spin,i) for i,count in enumerate(mu) for spin in (1,0)[:count]]
        expr=initial
        for step in steps:
            expr=apply(expr,derivative,step)
        check_highest(expr)
        weight,length=charge(expr)
        assert weight==(power-mu.count(2),power-sum(k>=1 for k in mu),*mu)
        assert length==2*power and not m.q(expr)
        polynomial=[[m.encode(tr),n] for tr,n in sorted(expr.items())]
        for step in reversed(steps):
            expr=apply(expr,reverse,step)
        assert list(expr)==[w]
        coefficient=orientation*expr[w]
        assert coefficient
        q_flavours=[weight[1]+k for k in mu]
        result.append(dict(weight=weight,length=length,level=6*power-len(steps),
                           flavours_for_interval=q_flavours,spin=weight[0]-weight[1],
                           terminal_grade=weight[1],steps=steps,reverse_coefficient=coefficient,
                           trace_terms=polynomial))
    assert len({tuple(v["weight"]) for v in result})==len(result)==10
    return sorted(result,key=lambda r:(-r["level"],r["weight"]))


def run():
    alphabet=tuple(l for l in letters_up_to(2,2) if 1<=l[0]+l[1]+len(l[2])<=2)
    assert len(alphabet)==17
    return dict(forward_classical_anticommutators=check_classical(alphabet),
                reverse_classical_anticommutators=check_reverse(alphabet),
                powers={str(p):records(p) for p in (5,6)})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    assert __debug__
    d=run()
    m.write(a.out,d)
    for p,rows in d['powers'].items():
        print(p,[(r['weight'],r['reverse_coefficient']) for r in rows])
