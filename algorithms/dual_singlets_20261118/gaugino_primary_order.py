"""Six ordered commutant actions on gaugino powers and the allowed primary source bidegrees."""
import argparse
from collections import Counter
from pathlib import Path

import mixed_incoming as m
from fortuity.full.letters import is_ferm
from rank_four_lowered_replay import derivative
from rank_four_sixth import G

STEPS=tuple((spin,flavour) for spin in range(2) for flavour in range(3))


def reverse(word,spin,flavour):
    result=Counter()
    phase=1
    for i,tr in enumerate(word):
        for j,letter in enumerate(tr):
            if letter==(0,0,(flavour,)):
                sg,out=m.so_canon(word[:i]+(tr[:j]+(G[spin],)+tr[j+1:],)+word[i+1:])
                if sg:
                    result[out]+=phase*sg
            if is_ferm(letter):
                phase=-phase
    return result


def apply(expr,action,step):
    result=Counter()
    for w,c in expr.items():
        for v,n in action(w,*step).items():
            result[v]+=n*c
    return {w:c for w,c in result.items() if c}


def charge(expr):
    values=set()
    for w in expr:
        ls=sum(w,())
        deg=[sum(x[i] for x in ls) for i in (0,1)]
        deg.extend(sum(i in x[2] for x in ls) for i in range(3))
        values.add((tuple(deg),len(ls)))
    assert len(values)==1
    return values.pop()


def run():
    result=[]
    for power in range(3,10):
        orientation,w=m.so_canon((G,)*power)
        expr={w:orientation}
        terms=[]
        for step in STEPS:
            expr=apply(expr,derivative,step)
            terms.append(len(expr))
        weight,length=charge(expr)
        assert weight==(power-3,power-3,2,2,2) and length==2*power
        assert not m.q(expr)
        primary=sorted([(m.encode(word),n) for word,n in expr.items()])
        for step in reversed(STEPS):
            expr=apply(expr,reverse,step)
        assert list(expr)==[w]
        coefficient=expr[w]*orientation
        assert coefficient!=0
        sources=[]
        for order in (1,2,3):
            lower=(weight[0]-order,weight[1]-order)+(weight[2]+order,)*3
            n=length+order-1
            if min(lower)>=0:
                assert sum(lower)-n==1 and m.sector(lower,n)==(1,1)+(2*(power-1),)*3
            sources.append(dict(order=order,derivative_flavour_degrees=lower,letters=n,permitted=min(lower)>=0))
        result.append(dict(power=power,steps=STEPS,terms_after_each_lowering=terms,primary_weight=weight,
                           primary_words=primary,reverse_coefficient=coefficient,primary_sources=sources))
    return result


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args();assert __debug__
    z=run();m.write(args.out,z)
    print([(r["power"],r["reverse_coefficient"]) for r in z])
