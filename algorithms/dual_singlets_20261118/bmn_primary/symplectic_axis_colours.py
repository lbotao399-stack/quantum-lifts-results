"""Keep exact sp(4) colour singlets on the first page of the pure-axis filtration."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from pathlib import Path

from symplectic_axis_bounds import partner, structure
from b2_schur_complex import multiplicity


@lru_cache(None)
def mixed_patterns(residue, length):
    """Count sets of distinct mixed Taylor indices by exterior-power multiplicities."""
    total = sum(residue)
    if length == 0:
        return {():1} if total==0 else {}
    if total<2*length or total-max(residue,default=0)<length:
        return {}
    jets = [w for w in product(*(range(n+1) for n in residue))
            if sum(v>0 for v in w)>=2 and sum(w)<=total-2*(length-1)]
    index = {w:i for i,w in enumerate(jets)}
    patterns = Counter()

    def visit(left, remaining, start, parts):
        if remaining==0:
            if not any(left):
                patterns[tuple(sorted(parts))] += 1
            return
        if sum(left)<2*remaining or sum(left)-max(left)<remaining:
            return
        if remaining==1:
            if left in index and index[left]>=start:
                patterns[tuple(sorted((*parts,1)))] += 1
            return
        for i in range(start,len(jets)):
            jet = jets[i]
            most = min([remaining]+[a//b for a,b in zip(left,jet) if b])
            for k in range(1,most+1):
                visit(tuple(a-k*b for a,b in zip(left,jet)), remaining-k, i+1, (*parts,k))

    visit(residue,length,0,())
    return patterns


def singlet_bound(axis_blocks):
    characters = defaultdict(Counter)
    degrees = defaultdict(set)
    for energy,(a,b),length,_,_,_,n in axis_blocks:
        if n:
            orbit = {(i*x,j*y) for x,y in permutations((a,b)) for i,j in product((-1,1),repeat=2)}
            for weight in orbit:
                characters["axis",energy,length][weight] = n
            degrees[energy].add(length)
    adjoint,_ = structure()
    for k in range(1,len(adjoint)+1):
        for selected in combinations(adjoint,k):
            weight = sum(x for x,y in selected),sum(y for x,y in selected)
            characters["mixed",k][weight] += 1

    @lru_cache(None)
    def tensor(factors):
        if not factors:
            return {(0,0):1}
        weights = Counter()
        for (a,b),n in tensor(factors[:-1]).items():
            for (x,y),k in characters[factors[-1]].items():
                weights[a+x,b+y] += n*k
        return weights

    @lru_cache(None)
    def singlets(factors):
        weights = tensor(factors)
        # Weyl alternation at rho=(2,1), in the defining symplectic weight basis.
        value = 0
        for (a,b),sign in (((2,1),1),((1,2),-1)):
            for i,j in product((-1,1),repeat=2):
                value += sign*i*j*weights.get((2-i*a,1-j*b),0)
        assert value>=0
        return value

    def at_weight(weight, length):
        count = 0
        choices = [[(e,y) for e in range(v+1) for y in degrees[e]] for v in weight]
        for axes in product(*choices):
            residue = tuple(sorted((v-e for v,(e,_) in zip(weight,axes) if v>e),reverse=True))
            mixed_length = length-sum(y for _,y in axes)
            if mixed_length<0 or sum(residue)<2*mixed_length:
                continue
            pure = [("axis",e,y) for e,y in axes if e]
            for pattern,n in mixed_patterns(residue,mixed_length).items():
                factors = tuple(sorted(pure+[("mixed",k) for k in pattern]))
                count += n*singlets(factors)
        return count
    return at_weight


def compute(axis_data):
    assert axis_data["group"]=="USp(4)" and axis_data["maximum_axis_energy"]==6
    bound = singlet_bound(axis_data["blocks"])
    assert bound((1,1),2)==1 and bound((1,1),1)==0
    rows = []
    for order,source_length,test_weight in ((2,13,(5,5,3,1)),(3,11,(6,6,1)),
                                            (3,13,(6,5,2)),(2,15,(4,4,4,2))):
        shape = partner((6,5,5),order)
        assert sum(test_weight)==sum(shape) and all(sum(test_weight[:k])<=sum(shape[:k])
                                                     for k in range(1,len(test_weight)+1))
        target_length = source_length+1-order
        full = multiplicity(shape,target_length)
        upper = bound(test_weight,target_length)
        rows.append(dict(order=order,source_length=source_length,target_diagram=shape,
                         target_length=target_length,test_weight=test_weight,target_cochains=full,
                         axis_page_singlets=upper,cohomology_upper=min(full,upper)))
    assert [r["cohomology_upper"] for r in rows]==[92,554,0,0]
    return dict(group="USp(4)",source_weight=[6,5,5],rows=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("axis_data",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    data = compute(json.loads(args.axis_data.read_bytes()))
    args.out.write_text(json.dumps(data,indent=2)+"\n")
