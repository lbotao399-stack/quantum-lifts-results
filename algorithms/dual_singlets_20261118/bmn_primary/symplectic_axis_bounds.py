"""Exact one-variable current homology of sp(4) and finite Taylor length bounds."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import sys

from flint import fmpz_mat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from b2_schur_complex import partner
from exterior_bound import partitions


def structure():
    basis, weights, coordinates = [], [], []
    units = ((1,0), (0,1))
    for i,j in product(range(2), repeat=2):
        matrix = fmpz_mat(4,4)
        matrix[i,j] = 1
        matrix[j+2,i+2] = -1
        basis.append(matrix)
        coordinates.append((i,j))
        weights.append(tuple(a-b for a,b in zip(units[i],units[j])))
    for upper in (True, False):
        for i,j in ((0,0),(0,1),(1,1)):
            matrix = fmpz_mat(4,4)
            x,y = ((0,2) if upper else (2,0))
            matrix[i+x,j+y] = matrix[j+x,i+y] = 1
            basis.append(matrix)
            coordinates.append((i+x,j+y))
            weights.append(tuple((a+b)*(1 if upper else -1) for a,b in zip(units[i],units[j])))
    terms = {}
    for i,j in product(range(10),repeat=2):
        bracket = basis[i]*basis[j]-basis[j]*basis[i]
        coeff = [int(bracket[u,v]) for u,v in coordinates]
        assert sum((c*m for c,m in zip(coeff,basis)), fmpz_mat(4,4))==bracket
        terms[i,j] = [(k,c) for k,c in enumerate(coeff) if c]
    return weights, terms


def weight_spaces(maximum, weights):
    by_energy = [[] for _ in range(maximum+1)]
    by_energy[0] = [((),0,0)]
    for power in range(1,maximum+1):
        for colour,(x,y) in enumerate(weights):
            letter = 10*(power-1)+colour
            for degree in range(maximum,power-1,-1):
                by_energy[degree].extend((word+(letter,),a+x,b+y)
                                         for word,a,b in by_energy[degree-power])
    result = defaultdict(list)
    for degree, states in enumerate(by_energy):
        for word,a,b in states:
            if a>=b>=0:
                result[degree,(a,b),len(word)].append(word)
    return result


def boundary(source, target, brackets):
    index = {word:i for i,word in enumerate(target)}
    matrix = fmpz_mat(len(target),len(source))
    for col, word in enumerate(source):
        for i,j in combinations(range(len(word)),2):
            left,right = word[i],word[j]
            neighbours = word[:i]+word[i+1:j]+word[j+1:]
            mode = left//10+right//10+1
            for colour,c in brackets[left%10,right%10]:
                inserted = 10*mode+colour
                if inserted in neighbours:
                    continue
                sign = (-1)**(i+j-1+sum(k<inserted for k in neighbours))
                row = index[tuple(sorted((*neighbours,inserted)))]
                matrix[row,col] += sign*c
    return matrix


def homology(maximum):
    weights, brackets = structure()
    spaces = weight_spaces(maximum,weights)
    records = []
    lengths = {d:set() for d in range(maximum+1)}
    total = Counter()
    for degree,weight in sorted({(d,w) for d,w,_ in spaces}):
        previous = fmpz_mat(0,0)
        previous_rank = 0
        for length in range(degree+2):
            source = spaces.get((degree,weight,length),())
            target = spaces.get((degree,weight,length-1),())
            matrix = boundary(source,target,brackets)
            assert (previous*matrix).is_zero()
            rank = matrix.rank()
            if target:
                dimension = len(target)-rank-previous_rank
                assert dimension >= 0
                records.append([degree,list(weight),length-1,len(target),rank,previous_rank,dimension])
                if dimension:
                    lengths[degree].add(length-1)
                    orbit_size = (1 if weight==(0,0) else 4 if 0 in weight or weight[0]==weight[1] else 8)
                    total[degree,length-1] += orbit_size*dimension
            previous,previous_rank = matrix,rank
    return records,lengths,[[d,y,n] for (d,y),n in sorted(total.items())]


def bounds(lengths):
    @lru_cache(None)
    def at_weight(weight):
        longest = -1
        for pure in product(*(tuple(d for d in range(w+1) if lengths[d]) for w in weight)):
            remaining = [w-d for w,d in zip(weight,pure)]
            n = sum(remaining)
            largest = max(remaining)
            if n and n==largest:
                continue
            mixed = min(n//2,n-largest)
            longest = max(longest,sum(max(lengths[d]) for d in pure)+mixed)
        return longest

    results = []
    for weight in ((6,5,4),(6,5,5),(7,4,4),(7,5,3),(7,5,4),(8,5,2)):
        for order in (2,3):
            shape = partner(weight,order)
            row = dict(source_weight=weight,order=order,target_diagram=shape)
            if shape is not None:
                # Kostka positivity guarantees that every dominated weight of
                # this size occurs in the Schur module on len(shape) even axes.
                candidates = (w for w in partitions(sum(shape),shape[0]) if len(w)<=len(shape)
                              and all(sum(w[:i])<=sum(shape[:i]) for i in range(1,len(w)+1)))
                cap, selected = min((at_weight(w),w) for w in candidates)
                row.update(test_weight=selected,maximum_target_length=cap,
                           maximum_source_length=cap+order-1)
            results.append(row)
    return results


def calculate():
    records,lengths,dimensions = homology(6)
    assert [sorted(lengths[i]) for i in range(7)] == [[0],[1],[2],[3],[3],[4],[]]
    return dict(group="USp(4)",maximum_axis_energy=6,
                block_columns=["energy","dominant_weight","length","dimension","incoming_rank",
                               "outgoing_rank","homology_dimension"],
                blocks=records,homology_dimensions=dimensions,taylor_bounds=bounds(lengths))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(calculate(),separators=(",",":"))+"\n")
