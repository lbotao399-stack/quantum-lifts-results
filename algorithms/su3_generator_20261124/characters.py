"""Independent exact SU(3) dimensions from symmetric-group characters.

The colour-conjugation character formula is the published B04 formula.
The word enumeration and the colour tensor character calculation share
neither polynomial evaluation matrices nor matrix-rank assumptions.
"""
from collections import Counter
from functools import cache
from math import comb,factorial
from pathlib import Path
import json,sys
from trace_algebra import contents,is_ferm

RUNTIME=Path(__file__).resolve().parents[1]/"su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/src"
sys.path.insert(0,str(RUNTIME))
from exact_dims import character,cycle_data,partitions


@cache
def endomorphism_traces(cycles,size):
    squared=tuple(sorted((m if m%2 else m//2 for m in cycles for _ in range(1 if m%2 else 2)),reverse=True))
    ordinary=transposed=0
    for shape in partitions(sum(cycles)):
        if len(shape)<=size:
            ordinary+=character(shape,cycles)**2
            transposed+=character(shape,squared)
    return ordinary,transposed


@cache
def adjoint_traces(cycles,size):
    choices={():1}
    for length,count in Counter(cycles).items():
        choices={old+(length,)*k:coefficient*(-1)**(count-k)*comb(count,k)
                 for old,coefficient in choices.items() for k in range(count+1)}
    ordinary=transposed=0
    for retained,coefficient in choices.items():
        full,reverse=endomorphism_traces(retained,size)
        ordinary+=coefficient*full;transposed+=coefficient*reverse
    return ordinary,(-1)**sum(cycles)*transposed


@cache
def odd_dimension(bosons,fermions,size):
    terms={():1};denominator=1
    for counts,odd in ((bosons,False),(fermions,True)):
        for count in counts:
            denominator*=factorial(count);following=Counter()
            for previous,coefficient in terms.items():
                for cycle,weight in cycle_data(count,odd):following[tuple(sorted(previous+cycle,reverse=True))]+=coefficient*weight
            terms={cycle:weight for cycle,weight in following.items() if weight}
    value=0
    for cycle,coefficient in terms.items():
        ordinary,conjugated=adjoint_traces(cycle,size);value+=coefficient*(ordinary-conjugated)
    assert value>=0 and value%(2*denominator)==0
    return value//(2*denominator)


def singlets(length,size=3):
    totals=[]
    for plus,minus in ((4,4),(5,3)):
        total=0
        for word in contents(length,plus,minus):
            counts=Counter(word)
            bosons=tuple(sorted((c for letter,c in counts.items() if not is_ferm(letter)),reverse=True))
            fermions=tuple(sorted((c for letter,c in counts.items() if is_ferm(letter)),reverse=True))
            total+=odd_dimension(bosons,fermions,size)
        totals.append(total)
    return dict(length=length,size=size,weight_dimensions=totals,dotted_singlets=totals[0]-totals[1])


if __name__=="__main__":
    result=[singlets(length,size) for length in (6,7,8) for size in (3,length)]
    assert [r["dotted_singlets"] for r in result]==[1169,2091,1896,4991,1588,6609]
    print(json.dumps(result,indent=2))
