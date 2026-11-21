"""Count SU(N) Casimir jet products when every letter carries the first flavour."""
from functools import lru_cache
from itertools import product


@lru_cache(None)
def weight_multiplicity(size, weight):
    if min(weight)<0:
        return 0
    plus,minus,length,second,third = weight
    zero=(0,0,0,0,0)
    multiplicities={zero:1}
    for n,i,j,b,c in product(range(2,min(length,size)+1),range(plus+1),
                             range(minus+1),range(min(1,second)+1),range(min(1,third)+1)):
        generator=(i,j,n,b,c)
        updated=multiplicities.copy()
        for original, coefficient in multiplicities.items():
            limit=min((m-k)//g for m,k,g in zip(weight,original,generator) if g)
            if (b+c)%2:
                limit=min(limit,1)
            for k in range(1,limit+1):
                target=tuple(x+k*g for x,g in zip(original,generator))
                updated[target]=updated.get(target,0)+coefficient
        multiplicities=updated
    return multiplicities.get(weight,0)


def highest_weight(size, weight, group):
    plus,minus,length,second,third=weight
    assert plus>=minus and length>=second>=third>=0 and group in ("SU+","SU-")
    if (group=="SU-") != (length%2==1):
        return 0
    # Weyl shifts that increase the first flavour are empty at length = first flavour.
    return sum((-1)**(spin+flavour)*weight_multiplicity(
        size,(plus+spin,minus-spin,length,second+flavour,third-flavour))
        for spin,flavour in product(range(2),repeat=2))
