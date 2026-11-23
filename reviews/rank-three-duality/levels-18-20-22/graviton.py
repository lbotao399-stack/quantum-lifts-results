"""Defining-trace Cartan restriction and divided-jet descendants from sp3_omega_products.py."""
from collections import Counter
from functools import lru_cache

def bits(mask):
    return tuple((mask>>j)&1 for j in range(3))


def parity(g):
    return (g[0]+g[1]+g[7].bit_count())%2


def weight(g):
    return (g[0]+g[5],g[1]+g[6],*(x+y for x,y in zip(g[2:5],bits(g[7]))))


def canonical(seq, odd):
    odd_items=[x for x in seq if odd(x)]
    if len(set(odd_items))!=len(odd_items):
        return 0,()
    return (-1)**sum(x>y for i,x in enumerate(odd_items) for y in odd_items[i+1:]),tuple(sorted(seq))


def letter_odd(l):
    return (1+l[2].bit_count())%2


@lru_cache(None)
def slot_odd(slot):
    return sum(letter_odd(l) for l in slot)%2


@lru_cache(maxsize=100000)
def combine_slots(first,second):
    return canonical(first+second,letter_odd)


@lru_cache(None)
def factor_polynomial(g):
    # On each diagonal the normalized graded symmetrization equals the ordered product.
    start=tuple((0,0,1<<i) for i,n in enumerate(g[2:5]) for _ in range(n))
    start+=((1,0,0),)*g[0]+((0,1,0),)*g[1]
    sign,slot=canonical(start,letter_odd)
    terms={slot:sign} if sign else {}
    outside=[2+j for j,n in enumerate(bits(g[7])) if n]+[0]*g[5]+[1]*g[6]
    for axis in outside:
        ans=Counter()
        for w,c in terms.items():
            phase=1
            for i,(a,b,theta) in enumerate(w):
                if axis<2:
                    l=(a+(axis==0),b+(axis==1),theta)
                    scale=(a,b)[axis]+1
                else:
                    bit=1<<(axis-2)
                    scale=0 if theta&bit else phase*(-1)**((theta//(2*bit)).bit_count())
                    l=(a,b,theta|bit)
                if scale:
                    sig,nw=canonical(w[:i]+(l,)+w[i+1:],letter_odd)
                    if sig:ans[nw]+=c*scale*sig
                if letter_odd((a,b,theta)):phase=-phase
        terms={w:c for w,c in ans.items() if c}
    return terms


def multiply_graviton(row, g):
    # Fold the labelled Cartan slots under simultaneous permutations.
    # As the multiplier is invariant this can be done after every factor.
    result=Counter()
    for slots,c in row.items():
        for word,n in factor_polynomial(g).items():
            for i in range(len(slots)):
                within,new=combine_slots(slots[i],word)
                if not within:continue
                crossing=(-1)**(slot_odd(word)*sum(slot_odd(v) for v in slots[i+1:]))
                sg,target=canonical(slots[:i]+(new,)+slots[i+1:],slot_odd)
                if sg:result[target]+=2*c*n*within*crossing*sg
    return {w:c for w,c in result.items() if c}


@lru_cache(maxsize=4500)
def cartan(row,rank=3):
    if not row:return {((),)*rank:1}
    return multiply_graviton(cartan(row[:-1],rank),row[-1])


