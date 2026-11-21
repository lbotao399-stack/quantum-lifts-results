"""Two-flavour cochain Euler indices and free Casimir product multiplicities."""
import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from b2_schur_complex import centralizer, partitions
from bmn_primary_index import colour
from exceptional_scalar_pairs import TARGETS
from exceptional_weyl_pairs import DEGREES


@lru_cache(None)
def multiplicity(group, x, y, fermions):
    entries = 2 if group == "E6" else 1
    values = [Fraction(0)]*entries
    if min(x,y) >= 0:
        for first,second,odd in product(partitions(x),partitions(y),partitions(fermions)):
            coefficient = Fraction((-1)**(fermions-len(odd)),
                                   centralizer(first)*centralizer(second)*centralizer(odd))
            for i,trace in enumerate(colour(tuple(sorted(first+second+odd)),group)):
                values[i] += coefficient*trace
    assert all(v.denominator==1 and v>=0 for v in values)
    return tuple(map(int,values))


def free_counts(group, a, b):
    # Store the ordinary weight multiplicities by length parity.
    gens = []
    for d in DEGREES[group]:
        if d>a+b:
            break
        gens.extend((d-k,k,d,False) for k in range(d+1))
        gens.extend((d-k,k+1,d,True) for k in range(d))
    @lru_cache(None)
    def count(i,x,y,parity):
        if not (x or y):
            return int(parity==0)
        if i == len(gens):
            return 0
        u,v,d,odd = gens[i]
        top = min([x+y]+([x//u] if u else [])+([y//v] if v else [])+([1] if odd else []))
        return sum(count(i+1,x-n*u,y-n*v,(parity-d*n)%2) for n in range(top+1))
    return [count(0,a,b,s)-count(0,a+1,b-1,s) for s in (0,1)]


def compute(group):
    a,b = TARGETS[group]
    rows = {}
    for k in range(b+1):
        direct = multiplicity(group,a-k,b-k,k)
        raised = multiplicity(group,a+1-k,b-1-k,k)
        counts = [u-v for u,v in zip(direct,raised)]
        assert min(counts)>=0
        rows[a+b-k] = counts
    euler = [sum((-1)**length*v[i] for length,v in rows.items()) for i in range(len(direct))]
    free = free_counts(group,a,b)
    assert euler[0]-free[0] == 1
    return {"group":group, "weight":[a,b,0], "lengths":rows,"ordinary_index":euler,
            "free_even_odd_length_multiplicities":free}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=[compute(g) for g in TARGETS]
    print(result)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
