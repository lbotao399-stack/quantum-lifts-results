"""Verify the rational product minors and two independent B2 characters at a BMN weight."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from math import lcm
from pathlib import Path

from flint import fmpz_mat, nmod_mat

import mixed_incoming as m
from sp3_omega_product_replay import check_factor, diagonal, odd_jet, odd_slot, ordered
from trace_descent import dimensions


def clean(terms):
    return {word:n for word,n in terms.items() if n}


@lru_cache(maxsize=300)
def labelled(factors):
    if not factors:
        return {((),()):1}
    result=Counter()
    for (left,right),n in labelled(factors[:-1]).items():
        for part,c in diagonal(factors[-1]).items():
            s,lr=ordered(left+part,odd_jet)
            if s:
                result[lr,right]+=2*n*c*s*(-1)**(odd_slot(part)*odd_slot(right))
            s,rr=ordered(right+part,odd_jet)
            if s:
                result[left,rr]+=2*n*c*s
    return clean(result)


@lru_cache(maxsize=36000)
def root_word(slots,row,column):
    result=Counter()
    for k,slot in enumerate(slots):
        for i,(_,_,mask) in enumerate(slot):
            if not (mask & (1<<column)) or mask & (1<<row):
                continue
            phase=(-1)**sum((mask>>j)&1 for j in range(min(row,column)+1,max(row,column)))
            jet=(0,0,mask^(1<<row)^(1<<column))
            sign,target=ordered(slot[:i]+(jet,)+slot[i+1:],odd_jet)
            if sign:
                result[slots[:k]+(target,)+slots[k+1:]]+=sign*phase
    return tuple(clean(result).items())


def action(terms,i,j):
    out=Counter()
    for w,n in terms.items():
        for v,c in root_word(w,i,j):
            out[v]+=c*n
    return clean(out)


def highest(terms,weight):
    # SU(3) extremal projector in the root order of the archived full_extremal.py.
    for i,j in ((1,2),(0,2),(0,1)):
        raised=terms
        result=Counter(terms)
        coefficient=Fraction(1)
        for power in range(1,sum(weight)+2):
            raised=action(raised,i,j)
            if not raised:
                break
            coefficient/=-power*(weight[i]-weight[j]+j-i+power)
            lowered=raised
            for _ in range(power):
                lowered=action(lowered,j,i)
            for w,n in lowered.items():
                result[w]+=coefficient*n
        else:
            raise AssertionError("Root string exceeds total polynomial degree")
        terms=clean(result)
    assert all(not action(terms,i,i+1) for i in (0,1))
    return terms


def witness(projected):
    denominators=[lcm(*(n.denominator for n in row.values())) for row in projected]
    assert all(v%m.P for v in denominators)
    columns=sorted({w for row in projected for w in row})
    integral=[[int(d*row.get(w,0)) for w in columns] for row,d in zip(projected,denominators)]
    pivots=m.pivotcols(nmod_mat(integral,m.P))
    assert len(pivots)==len(integral)
    square=fmpz_mat([[row[p] for p in pivots] for row in integral])
    det=int(square.det())
    assert det%m.P
    return dict(rank=len(integral),denominators=denominators,signed_integer_determinant=str(det),
                selected_labelled_columns=[columns[i] for i in pivots])


def run(path, weight=(10,4,4)):
    raw=path.read_bytes()
    data=json.loads(gzip.decompress(raw) if path.suffix==".gz" else raw)
    weight=tuple(weight)
    if isinstance(data,list):
        data,=[entry for entry in data if tuple(entry["weight"])==weight]
    assert (data["group"],data["level"],tuple(data["weight"]))==("USp(4)",2*sum(weight),weight)
    char={}
    for n in range(max(weight),sum(weight)+1):
        # The orthogonal colour formula uses the exterior square of the defining
        # five-dimensional representation. The symplectic formula uses Sym² of four.
        orthogonal=dimensions((0,0,*weight),n,5,"SO")
        symplectic=dimensions((0,0,*weight),n,4,"Sp")
        assert orthogonal==symplectic and sum(orthogonal.values())==data["dimensions"][str(n)]
        char[n]=orthogonal
    minors={}
    for record in data["even_length_records"]:
        products=[tuple(map(tuple,w)) for w in record["products"]]
        if products:
            for g in {g for w in products for g in w}:
                check_factor(g)
            rows=[highest(labelled(w),weight) for w in products]
            minors[record["length"]]=witness(rows)
            assert len(rows)==record["product_lower_bound"]
            print(record["length"],len(rows),minors[record["length"]]["signed_integer_determinant"],flush=True)
    even=sum(sum(vals.values()) for n,vals in char.items() if n%2==0)
    odd=sum(sum(vals.values()) for n,vals in char.items() if n%2)
    known=sum(d["rank"] for d in minors.values())
    assert even-odd==data["euler_with_length_sign"] and known==data["even_product_lower_bound"]
    bound=max(0,known-even+odd)
    assert bound==data["odd_cohomology_lower_bound"]
    return dict(index_data_sha256=hashlib.sha256(raw).hexdigest(),weight=weight,
                cochains_by_length_and_fermion_degree=char,even_cochains=even,odd_cochains=odd,
                even_product_minors=minors,odd_cohomology_lower_bound=bound)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,required=True)
    parser.add_argument("--weight",type=int,nargs=3,default=(10,4,4))
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    m.write(args.out,run(args.data,args.weight))
