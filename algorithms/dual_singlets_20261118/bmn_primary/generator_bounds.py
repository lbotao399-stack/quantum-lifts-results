"""Test whether specified odd BMN generators can supply the indexed primary spaces."""
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from pathlib import Path


@lru_cache(None)
def character(shape, with_descendants):
    a,b,c=shape
    weights=Counter()
    for i in range(b,a+1):
        for j in range(c,b+1):
            for k in range(j,i+1):
                weights[k,i+j-k,a+b+c-i-j]+=1
    if not with_descendants:
        return weights
    total=Counter()
    for w,n in weights.items():
        for extra in product(range(2),repeat=3):
            total[tuple(x+y for x,y in zip(w,extra))]+=n
    return total


@lru_cache(maxsize=20000)
def dressings(target, core, with_descendants=False):
    """Return the nonnegative Littlewood--Richardson coefficients of all dressings."""
    if any(i<j for i,j in zip(target,core)):
        return {}
    shifted=target[0]+2,target[1]+1,target[2]
    result=Counter()
    for (a,b,c),n in character(core,with_descendants).items():
        x,y,z=shifted[0]-a,shifted[1]-b,shifted[2]-c
        if min(x,y,z)<0 or x==y or x==z or y==z:
            continue
        sign=(-1)**((x<y)+(x<z)+(y<z))
        i,j,k=sorted((x,y,z),reverse=True)
        result[i-2,j-1,k]+=sign*n
    if any(n<0 for n in result.values()):
        raise ArithmeticError("Negative tensor-product multiplicity")
    return {w:n for w,n in result.items() if n}


def table(data):
    return {tuple(row[:3]):row[3] for row in data["rows"]}


def first_failure(cores, ordinary, primary):
    """Bound the even dressing recursively under the proposed list of odd generators."""
    @lru_cache(None)
    def even(weight):
        count=ordinary[weight]
        for core,m in cores.items():
            count+=m*sum(n*even(w) for w,n in dressings(weight,core,True).items())
        return count

    for target,index in sorted(primary.items(),key=lambda item:(sum(item[0]),item[0])):
        if index>=0:
            continue
        parts=[]
        for core,m in cores.items():
            terms=dressings(target,core)
            parts.append(m*sum(n*even(w) for w,n in terms.items()))
        if sum(parts)<-index:
            return dict(target=target,level=2*sum(target),required=-index,upper_by_core=parts)
    return None


def compare(ordinary_data, primary_data):
    group=ordinary_data["group"]
    if primary_data["group"]!=group or ordinary_data["index_type"]!="ordinary" or primary_data["index_type"]!="primary":
        raise ValueError("Ordinary and primary tables do not match")
    ordinary,primary=table(ordinary_data),table(primary_data)
    if group=="G2":
        cores=[(6,5,4),(6,6,4),(6,6,6)]
        baseline={}
    elif group=="C2":
        baseline={(10,4,4):2}
        # All dressings entering these six tests have degree at most 14, where
        # the 297 adjacent-rank equalities establish odd vanishing over Q.
        early={w:q for w,q in primary.items() if sum(w)<=29}
        cores=[]
        for a in range(9):
            for b in range(a+1):
                for c in range(b+1):
                    core=a,b,c
                    if sum(core)>=15 and first_failure({**baseline,core:1},ordinary,early) is None:
                        cores.append(core)
    else:
        raise ValueError("This comparison uses the G2 and C2 preceding-rank certificates")
    trials=[]
    for core in cores:
        proposed={**baseline,core:1}
        trials.append(dict(core=core,first_failure=first_failure(proposed,ordinary,primary)))
    return dict(group=group,baseline=[[w,n] for w,n in baseline.items()],trials=trials)


if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("ordinary",type=Path)
    p.add_argument("primary",type=Path)
    args=p.parse_args()
    print(json.dumps(compare(json.loads(args.ordinary.read_bytes()),json.loads(args.primary.read_bytes())),indent=2))
