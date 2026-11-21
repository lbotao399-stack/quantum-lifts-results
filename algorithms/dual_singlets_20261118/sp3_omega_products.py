"""Project explicit single-graviton products in the prospective incoming Sp(3) source."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement, product
from pathlib import Path

import numpy as np
import mixed_incoming as m
from mixed33_exact import apply, exact
from orthosymplectic_classical_page import minor


WEIGHT=(1,1,7,2,2)
LENGTH=10


def bits(mask):
    return tuple((mask>>j)&1 for j in range(3))


def parity(g):
    return (g[0]+g[1]+g[7].bit_count())%2


def family(g):
    return (sum(g[:5]),sum(g[:2]),sum(g[5:7]),g[7].bit_count())


def weight(g):
    return (g[0]+g[5],g[1]+g[6],*(x+y for x,y in zip(g[2:5],bits(g[7]))))


def canonical(seq, odd):
    odd_items=[x for x in seq if odd(x)]
    if len(set(odd_items))!=len(odd_items):
        return 0,()
    return (-1)**sum(x>y for i,x in enumerate(odd_items) for y in odd_items[i+1:]),tuple(sorted(seq))


def all_products(max_power, target_weight=WEIGHT, length=LENGTH, *, include_odd_powers=False):
    by_key=defaultdict(list)
    for k in range(2,max_power+1,1 if include_odd_powers else 2):
        for hp,hm,sp,sm,mask in product(range(2),range(2),range(target_weight[0]+1),
                                        range(target_weight[1]+1),range(8)):
            for a in range(k-hp-hm+1):
                for b in range(k-hp-hm-a+1):
                    f=(hp,hm,a,b,k-hp-hm-a-b,sp,sm,mask)
                    w=weight(f)
                    if all(x<=y for x,y in zip(w,target_weight)) and sum(w)-k<=sum(target_weight)-length:
                        by_key[(k,w)].append(f)
    keys=sorted(by_key)
    def types(length,remaining,start):
        if length==0:
            if not any(remaining):
                yield ()
            return
        for i in range(start,len(keys)):
            k,w=keys[i]
            if k>length:break
            rest=tuple(y-x for x,y in zip(w,remaining))
            if min(rest)<0 or sum(rest)<length-k:continue
            for tail in types(length-k,rest,i):
                yield (keys[i],)+tail
    grouped=defaultdict(list)
    for pattern in types(length,target_weight,0):
        repeats=Counter(pattern)
        choices=[combinations_with_replacement(by_key[key],n) for key,n in repeats.items()]
        for terms in product(*choices):
            row=tuple(sorted(g for part in terms for g in part))
            if len({g for g in row if parity(g)})!=sum(parity(g) for g in row):continue
            grouped[tuple(sorted(family(g) for g in row))].append(row)
    return grouped


def raised_factor(g,axis):
    # Simple roots of dotted spin, then the two adjacent flavour roots.
    if axis==0:
        if g[1] and not g[0]:
            yield 1,(1,0,*g[2:])
        if g[6]:
            yield g[6],(*g[:5],g[5]+1,g[6]-1,g[7])
    else:
        i=axis-1
        if g[3+i]:
            v=list(g)
            v[2+i]+=1
            v[3+i]-=1
            yield g[3+i],tuple(v)
        if g[7]&(2<<i) and not g[7]&(1<<i):
            yield 1,(*g[:7],g[7]^(3<<i))


def raising(row):
    out=Counter()
    for i,g in enumerate(row):
        for axis in range(3):
            for n,h in raised_factor(g,axis):
                sign,word=canonical(row[:i]+(h,)+row[i+1:],parity)
                if sign:out[axis,word]+=n*sign
    return {v:c for v,c in out.items() if c}


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


def run(max_power,out,stop,target_weight=WEIGHT,length=LENGTH,rank=3):
    families=all_products(max_power,target_weight,length)
    print("raw products",sum(map(len,families.values())),"families",len(families),
          "largest",max(map(len,families.values())),flush=True)
    basis,recipes=[],[]
    echelon={}
    for key in sorted(families,key=lambda k:(max(i[0] for i in k),k)):
        words=families[key]
        images=[raising(w) for w in words]
        _,_,kernel=exact(images)
        for combination in kernel:
            assert not apply(combination,images)
            expr=Counter()
            for i,c in combination.items():
                for w,n in cartan(words[i],rank).items():expr[w]+=c*n
            integral={w:n for w,n in expr.items() if n}
            reduced={w:n%m.P for w,n in integral.items() if n%m.P}
            while reduced:
                lead=min(reduced)
                if lead not in echelon:
                    inv=pow(reduced[lead],-1,m.P)
                    echelon[lead]={w:n*inv%m.P for w,n in reduced.items()}
                    basis.append(integral)
                    recipes.append([[words[i],n] for i,n in sorted(combination.items())])
                    break
                mul=reduced[lead]
                for w,n in echelon[lead].items():
                    t=(reduced.get(w,0)-mul*n)%m.P
                    if t:reduced[w]=t
                    else:reduced.pop(w,None)
        if stop and len(basis)>=stop:break
    columns=sorted(echelon)
    matrix=np.asarray([[row.get(k,0)%m.P for k in columns] for row in basis],dtype=np.int64)
    witness=minor(matrix,m.P)
    assert witness["rank"]==len(basis)
    result=dict(weight=target_weight,length=length,group=f"USp({2*rank})",defining_trace=True,
                maximum_trace_power=max_power,prime=m.P,cartan_rank=len(basis),
                minor=witness,cartan_columns=columns,highest_weight_recipes=recipes)
    m.write(out,result)
    print("protected product minor",len(basis),witness["determinant"],flush=True)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-power",type=int,default=6,choices=(2,4,6,8,10))
    parser.add_argument("--stop",type=int,default=454)
    parser.add_argument("--weight",type=int,nargs=5,default=WEIGHT)
    parser.add_argument("--length",type=int,default=LENGTH)
    parser.add_argument("--rank",type=int,default=3)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    run(args.max_power,args.out,args.stop,tuple(args.weight),args.length,args.rank)
