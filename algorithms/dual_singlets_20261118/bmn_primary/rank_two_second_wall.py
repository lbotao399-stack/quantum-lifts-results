"""Bound Cartan kernels of the USp(4) second target on its proper connected subgroups."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from charge_maps import content_patterns
from generator_bounds import dressings
from su3_levi_branching import branches
from su3_levi_kernel import check, key, resolve
from su3_third_levi_weights import read_records
from trace_descent import dimensions


@lru_cache(None)
def active_dimension(weight, length):
    if length == 0:
        return int(not any(weight))
    if not max(weight[:3]) <= length <= sum(weight):
        return 0
    return sum(dimensions((*weight[3:], *weight[:3]), length, 2, "Sp").values())


def tensor_blocks(weight, length):
    """Unordered SU(2) factors. Return multiplicity before symmetrizing in colour."""
    for first in product(*(range(k+1) for k in weight[:3])):
        if not first[0] >= first[1] >= first[2]:
            continue
        for second, copies in dressings(weight[:3], first).items():
            # (i,j) times (k,l) contains (i+k-h,j+l+h), 0 <= h <= min(i-j,k-l).
            for i in range(weight[3]+1):
                for j in range(min(i,weight[4])+1):
                    for h in range(i-j+1):
                        k, l = weight[3]-i+h, weight[4]-j-h
                        if l<0 or k-l<h:
                            continue
                        a, b = first+(i,j), second+(k,l)
                        for n in range(length+1):
                            if ((a,n) <= (b,length-n) and active_dimension(a,n)
                                    and active_dimension(b,length-n)):
                                yield a,n,b,length-n,copies


@lru_cache(None)
def power(k, exterior):
    # Weights in units of the roots of the two commuting sl(2) factors.
    weights = ((-1,0),(0,0),(1,0),(0,-1),(0,0),(0,1))
    choose = combinations if exterior else combinations_with_replacement
    return Counter((sum(x for x,y in row),sum(y for x,y in row)) for row in choose(weights,k))


@lru_cache(None)
def colour_dimension(bosons, fermions):
    weights = {(0,0):1}
    for p, exterior in [(k,False) for k in bosons]+[(k,True) for k in fermions]:
        result = Counter()
        for (x,y),n in power(p,exterior).items():
            for (a,b),m in weights.items():
                result[a+x,b+y] += m*n
        weights = result
    value = weights.get((0,0),0)-weights.get((1,0),0)-weights.get((0,1),0)+weights.get((1,1),0)
    assert value>=0
    return value


def direct_cochains(weight, length):
    """Independent Weyl integration on the six-dimensional direct-sum adjoint."""
    count = 0
    for ds in (0,1):
        for perm in permutations(range(3)):
            shifted = (weight[3]+ds,weight[4]-ds) + tuple(weight[j]+perm[j]-j for j in range(3))
            sign = (-1)**(ds + sum(i>j for k,i in enumerate(perm) for j in perm[k+1:]))
            for bosons,fermions,n in content_patterns(shifted,length):
                count += sign*n*colour_dimension(tuple(bosons),tuple(fermions))
    return count


def calculate(directory):
    weight, length, prime = (4,3,3,2,2),10,2147483647
    old = HERE.parent/"taylor_hooks/su3_scalar_levi/full_block"
    paths = [old/name for name in ("minors_2147483647.json.gz","third_remaining_minors.json.gz",
                                    "third_single_flavour_minors.json.gz")]
    local = directory/"su2_minors.json.gz"
    receipt = json.loads((directory/"receipt.json").read_bytes())
    assert receipt == dict(job="7658782",tasks=65,files={local.name:hashlib.sha256(local.read_bytes()).hexdigest()})
    known = read_records(paths)
    fresh = read_records([local])
    assert not (set(fresh) & set(known)) and len(fresh)==65
    records = {**known,**fresh}
    used = set()

    def counts(w,n):
        if not n:
            assert not any(w)
            return 1,1,0
        label = key(w,n)
        used.add(label)
        value = check(records[label],w,n,prime)
        assert active_dimension(w,n) == value[0]
        return value

    plan = branches(weight,length,1)
    for row in plan:
        counts(row["weight"],row["length"])
    wall = resolve({"wall":dict(weight=weight,length=length,colour_parity=1,active=plan)},records,prime)["sectors"]["wall"]
    totals = [0,0,0]
    blocks = list(tensor_blocks(weight,length))
    for a,i,b,j,m in blocks:
        da,ha,ka = counts(a,i)
        db,hb,kb = counts(b,j)
        ordered = m*(1 if (a,i)==(b,j) else 2)
        for k, v in enumerate((da*db, ha*hb, ka*hb+ha*kb-ka*kb)):
            totals[k] += ordered*v
    principal = counts(weight,length)
    assert set(fresh)<=used and len(used)==529
    assert totals[0] == direct_cochains(weight,length)
    assert (totals,wall["cochains"],wall["classical_dimension_upper"],wall["cartan_kernel_upper"]) == (
        [137907,694,0],24078,726,0)
    assert principal[2]==0
    return dict(group="USp(4)",source_weight=[6,5,5],source_length=11,order=2,
                target_diagram=[5,5,3,1],containing_weight=weight,target_length=length,prime=prime,
                active_sectors=len(used),previous_records=len(used)-len(fresh),new_records=len(fresh),
                levi=wall,product=dict(unordered_blocks=len(blocks),cochains=totals[0],
                                      classical_dimension_upper=totals[1],cartan_kernel_upper=totals[2]),
                principal=dict(cochains=principal[0],classical_dimension_upper=principal[1],cartan_kernel_upper=principal[2]))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:
        parser.error("enable the rank and character checks")
    record=calculate(HERE/"7658782")
    args.out.write_text(json.dumps(record,indent=2)+"\n")
