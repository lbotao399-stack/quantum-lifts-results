"""Count the remaining scalar-commutator equations after exact ghost commutations in traces."""
from collections import Counter
from functools import cache
from itertools import permutations, product
import json
from math import comb


def sign(perm):
    return (-1)**sum(x>y for i,x in enumerate(perm) for y in perm[i+1:])


def rotate(word):
    return min((word[i:]+word[:i] for i in range(len(word))),default=())


def generators(scalars_commute):
    out = Counter()
    for length in range(7):
        for word in product(range(3),repeat=length):
            weight = tuple(word.count(j) for j in range(3))
            for u,v in product(range(2),repeat=2):
                if u+v+length < 2:
                    continue
                if scalars_commute:
                    if word != tuple(sorted(word)) or (length+u+v)%2:
                        continue
                else:
                    forward = rotate(word)
                    backward = rotate(word[::-1])
                    if word != min(forward,backward):
                        continue
                    if forward == backward and (length+u+v)%2:
                        continue
                out[u,v,*weight] += 1
    return out


@cache
def trace_counts(scalars_commute):
    dp = Counter({(0,0,0,0,0,0):1})
    for weight, multiplicity in sorted(generators(scalars_commute).items()):
        odd = sum(weight[:2])%2
        updated = dp.copy()
        for key,n in dp.items():
            max_count = min([(cap-x)//y for x,y,cap in zip(key,weight,(5,4,4,4,4)) if y] + [7-key[-1]])
            if sum(weight[2:]):
                max_count = min(max_count,(6-sum(key[2:5]))//sum(weight[2:]))
            if odd:
                max_count = min(max_count,multiplicity)
            for k in range(1,max_count+1):
                target = tuple(x+k*y for x,y in zip(key,weight))+(key[-1]+k,)
                coefficient = comb(multiplicity,k) if odd else comb(multiplicity+k-1,k)
                updated[target] += n*coefficient
        dp = updated
    return dp


def highest(counts,dotted=4):
    result = Counter()
    for spin in (0,1):
        for perm in permutations(range(3)):
            scalar = tuple(2+perm[i]-i for i in range(3))
            weight = (dotted+spin,dotted-spin,*scalar)
            for traces in range(1,8):
                result[traces] += (-1)**spin*sign(perm)*counts[(*weight,traces)]
    assert all(n>=0 for n in result.values())
    return result


def calculate():
    partial,full = [highest(trace_counts(mode)) for mode in (False,True)]
    rows = [dict(traces=k,ghost_relations=partial[k],all_relations=full[k],
                 remaining_scalar_commutator_dimension=partial[k]-full[k]) for k in range(1,8)]
    assert all(row["remaining_scalar_commutator_dimension"] >= 0 for row in rows)
    divisible = []
    for power in range(5):
        p, f = [highest(trace_counts(mode),4-power) for mode in (False,True)]
        divisible.append(sum(p.values())-sum(f.values()))
    exact = [x-y for x,y in zip(divisible,divisible[1:]+[0])]
    return dict(dotted_degrees=[4,4],flavour_degrees=[2,2,2],letter_length=14,
                by_trace_number=rows,partial_dimension=sum(partial.values()),
                full_dimension=sum(full.values()),scalar_commutator_dimension=sum(partial.values())-sum(full.values()),
                dimensions_divisible_by_glueball_power=divisible,
                dimensions_by_exact_glueball_power=exact)


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
