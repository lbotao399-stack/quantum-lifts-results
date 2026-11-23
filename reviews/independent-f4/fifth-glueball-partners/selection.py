"""The physical singlet copy of (4,4,3) and the three possible Ward partners."""
from itertools import permutations
import json


def sign(order):
    return (-1)**sum(a>b for i,a in enumerate(order) for b in order[i+1:])


def product(left,right):
    answer={}
    for a,x in left.items():
        for b,y in right.items():
            if a&b:
                continue
            crosses=sum((b & ((1<<i)-1)).bit_count() for i in range(16) if a&(1<<i))
            mask=a|b
            answer[mask]=answer.get(mask,0)+(-1)**crosses*x*y
    return {k:v for k,v in answer.items() if v}


def determinant(n,scalar):
    result={}
    for order in permutations(range(n)):
        term={0:sign(order)}
        for i,j in enumerate(order):
            entry={(1<<i)|(1<<(8+j)):1,(1<<(4+i))|(1<<(12+j)):1}
            if i==j and i<3:
                entry[0]=scalar[i]
            term=product(term,entry)
        for mask,coeff in term.items():
            result[mask]=result.get(mask,0)+coeff
    return {k:v for k,v in result.items() if v}


def run():
    coefficients=[]
    for scalar in ((1,1,1),(2,3,5)):
        four=determinant(4,scalar)
        value=product(product(four,four),determinant(3,scalar))[(1<<16)-1]
        coefficients.append(value)
    assert coefficients==[14400,432000]
    rows=[]
    for t in range(6):
        odd=5-t
        if odd and t<3:
            continue
        shape=(t,t)+(3,)*odd
        assert sum(shape)==15-t
        rows.append(dict(dotted_each=t,flavour_each=odd,letters=14-t,diagram=shape))
    assert [tuple(r["diagram"]) for r in rows]==[(3,3,3,3),(4,4,3),(5,5)]
    return dict(polarization_coefficients=coefficients,physical_partner_diagrams=rows)


if __name__=="__main__":
    print(json.dumps(run(),indent=2))
