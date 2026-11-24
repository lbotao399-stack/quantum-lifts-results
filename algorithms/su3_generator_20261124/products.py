"""Actual closed trace polynomials for every graviton recipe in this grading."""
from collections import Counter,defaultdict
from functools import cache
from fractions import Fraction
from itertools import product,permutations,combinations_with_replacement
from math import factorial,gcd
from pathlib import Path
import argparse,json,pickle
from flint import fmpz_mat
from trace_algebra import canonical_multi,is_ferm,differential,raise_row,conjugate


def canonical(items,odd):
    odd_items=[item for item in items if odd(item)]
    if len(set(odd_items))!=len(odd_items):return 0,()
    return (-1)**sum(x>y for i,x in enumerate(odd_items) for y in odd_items[i+1:]),tuple(sorted(items))


def parity(g):return (g[0]+g[1]+g[7].bit_count())%2
def family(g):return (sum(g[:5]),sum(g[:2]),sum(g[5:7]),g[7].bit_count())
def weight(g):return (g[0]+g[5],g[1]+g[6],*(g[i+2]+((g[7]>>i)&1) for i in range(3)))


def all_products(maximum=7,plus=4,minus=4):
    target=(plus,minus,3,0,0);by_key=defaultdict(list)
    for length in range(2,maximum+1):
        for hp,hm,sp,sm,mask in product(range(2),range(2),range(plus+1),range(minus+1),range(2)):
            spec=(hp,hm,length-hp-hm,0,0,sp,sm,mask)
            if spec[2]<0:continue
            charges=weight(spec)
            if all(a<=b for a,b in zip(charges,target)) and sum(charges)-length<=4:
                by_key[length,charges].append(spec)
    keys=sorted(by_key)
    def choose(length,remaining,start):
        if not length:
            if not any(remaining):yield ()
            return
        for i in range(start,len(keys)):
            degree,charges=keys[i]
            if degree>length:break
            rest=tuple(b-a for a,b in zip(charges,remaining))
            if min(rest)<0 or sum(rest)<length-degree:continue
            for tail in choose(length-degree,rest,i):yield (keys[i],)+tail
    grouped=defaultdict(list)
    for selected in choose(7,target,0):
        choices=[combinations_with_replacement(by_key[key],count) for key,count in Counter(selected).items()]
        for choices_by_type in product(*choices):
            recipe=tuple(sorted(g for values in choices_by_type for g in values))
            if len({g for g in recipe if parity(g)})!=sum(map(parity,recipe)):continue
            grouped[tuple(sorted(map(family,recipe)))].append(recipe)
    return grouped


def raised_recipe(recipe):
    result=Counter()
    for i,g in enumerate(recipe):
        changed=[]
        if g[1] and not g[0]:changed.append((1,(1,0,*g[2:])))
        if g[6]:changed.append((g[6],(*g[:5],g[5]+1,g[6]-1,g[7])))
        for coefficient,value in changed:
            sign,target=canonical(recipe[:i]+(value,)+recipe[i+1:],parity)
            if sign:result[target]+=coefficient*sign
    return {k:v for k,v in result.items() if v}


def derivation(row,axis):
    result=Counter()
    for multi,c in row.items():
        phase=1
        for i,trace in enumerate(multi):
            for j,(a,b,theta) in enumerate(trace):
                if axis<2:
                    coefficient=(a,b)[axis]+1
                    letter=(a+(axis==0),b+(axis==1),theta)
                else:
                    bit=axis-2
                    coefficient=0 if bit in theta else phase*(-1)**sum(t>bit for t in theta)
                    letter=(a,b,tuple(sorted(theta+(bit,))))
                if coefficient:
                    changed=trace[:j]+(letter,)+trace[j+1:]
                    sign,target=canonical_multi(multi[:i]+(changed,)+multi[i+1:])
                    if sign:result[target]+=c*coefficient*sign
                if is_ferm((a,b,theta)):phase=-phase
    return {w:c for w,c in result.items() if c}


@cache
def factor(spec):
    letters=tuple((0,0,(i,)) for i,count in enumerate(spec[2:5]) for _ in range(count))
    letters+=((1,0,()),)*spec[0]+((0,1,()),)*spec[1]
    result=Counter();length=len(letters)
    for order in permutations(range(length)):
        odd=[j for j in order if is_ferm(letters[j])]
        phase=(-1)**sum(a>b for i,a in enumerate(odd) for b in odd[i+1:])
        sign,word=canonical_multi((tuple(letters[j] for j in order),))
        if sign:result[word]+=Fraction(phase*sign,factorial(length))
    result={w:c for w,c in result.items() if c}
    for axis in [2+i for i in range(3) if spec[7]&(1<<i)]+[0]*spec[5]+[1]*spec[6]:result=derivation(result,axis)
    assert not differential(result)
    return result


@cache
def recipe_polynomial(recipe):
    if not recipe:return {():Fraction(1)}
    result=Counter()
    for left,c in recipe_polynomial(recipe[:-1]).items():
        for right,value in factor(recipe[-1]).items():
            sign,word=canonical_multi(left+right)
            if sign:result[word]+=c*value*sign
    return {w:c for w,c in result.items() if c}


def combination(terms):
    result=Counter()
    for recipe,c in terms:
        recipe=tuple(tuple(g) for g in recipe)
        for word,value in recipe_polynomial(recipe).items():result[word]+=c*value
    row={w:c for w,c in result.items() if c}
    assert not differential(row) and not raise_row(row)
    assert conjugate(row)=={w:-c for w,c in row.items()}
    return row


def complete_gravitons():
    groups=all_products();assert sum(map(len,groups.values()))==447
    assert sum(map(len,all_products(7,5,3).values()))==361
    answer=[];recipes_out=[]
    for key,recipes in sorted(groups.items()):
        images=list(map(raised_recipe,recipes));keys=sorted({k for r in images for k in r});ids={k:i for i,k in enumerate(keys)}
        matrix=fmpz_mat(len(keys),len(images))
        for j,row in enumerate(images):
            for word,c in row.items():matrix[ids[word],j]=c
        kernel,dimension=matrix.nullspace()
        for j in range(dimension):
            coefficients={i:int(kernel[i,j]) for i in range(len(recipes)) if kernel[i,j]}
            divisor=gcd(*coefficients.values());terms=[(recipes[i],c//divisor) for i,c in coefficients.items()]
            answer.append(combination(terms));recipes_out.append(terms)
    return answer,recipes_out


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("archived",type=Path);parser.add_argument("output",type=Path)
    args=parser.parse_args();archive=json.loads(args.archived.read_text())
    proper=[combination(terms) for terms in archive["highest_weight_recipes"]]
    complete,recipes=complete_gravitons()
    output=dict(proper=proper,complete=complete,complete_recipes=recipes)
    args.output.write_bytes(pickle.dumps(output,protocol=5))
    print(json.dumps(dict(proper=len(proper),complete=len(complete),proper_terms=sum(map(len,proper)),complete_terms=sum(map(len,complete)))),flush=True)
