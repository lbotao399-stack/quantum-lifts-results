"""Project a two-variable positive-current complex onto the Taylor determinant before taking ranks."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import factorial
from exact_rank import rational_rank

@cache
def determinant_packet(shape):
    """Eliminate the image of y d/dx within a fixed colour and jet-degree occupation."""
    total = sum(degree * count for degree, count in shape)
    assert total % 2 == 0
    middle = total // 2
    scale = factorial(total)
    choices = [tuple(combinations(range(degree+1), count)) for degree, count in shape]
    states = defaultdict(list)
    for groups in product(*choices):
        state = sum(groups, ())
        if sum(state) in (middle, middle+1):
            states[sum(state)].append(state)
    rows = {}
    for state in states[middle+1]:
        row = {}
        start = 0
        for degree, count in shape:
            for k in range(start, start+count):
                a = state[k]
                if a and (k == start or a-1 != state[k-1]):
                    target = state[:k]+(a-1,)+state[k+1:]
                    row[target] = Fraction(a)
            start += count
        while row:
            lead = min(row)
            if lead not in rows:
                unit = row[lead]
                rows[lead] = {k:v/unit for k,v in row.items() if k != lead}
                break
            n = row.pop(lead)
            for k, c in rows[lead].items():
                delta = row.get(k, 0) - n*c
                if delta:
                    row[k] = delta
                else:
                    row.pop(k, None)
        else:
            raise ArithmeticError('Taylor lowering lost a positive-weight column')
    free = tuple(state for state in states[middle] if state not in rows)
    assert len(free) == len(states[middle])-len(states[middle+1])
    ids = {state:k for k,state in enumerate(free)}
    projected = {}
    for state in sorted(states[middle], reverse=True):
        if state in ids:
            projected[state] = {ids[state]:scale}
            continue
        value = Counter()
        for k,n in rows[state].items():
            for j,c in projected[k].items():
                value[j] -= n*c
        assert all(n.denominator == 1 for n in value.values())
        projected[state] = {j:int(n) for j,n in value.items() if n}
    return free, projected

def packet(letters):
    counts = Counter((d,c) for d,c,a in letters)
    groups = tuple(sorted((d,n,c) for (d,c),n in counts.items()))
    order = [(d,counts[d,c],c,a) for d,c,a in letters]
    sign = (-1)**sum(a>b for i,a in enumerate(order) for b in order[i+1:])
    state = tuple(a for d,n,c,a in sorted(order))
    return groups,state,sign

def basis(ce, length, weight):
    """Enumerate central-weight words once to find every nonempty occupation."""
    letters = [(a+b,c,a) for a,b in ce.modes for c in range(len(ce.weights))]
    found=set()
    for word in ce.words(length,weight):
        counts=Counter((letters[i][0],letters[i][1]) for i in word)
        found.add(tuple(sorted((d,n,c) for (d,c),n in counts.items())))
    words, projection = [], {}
    for groups in sorted(found):
        shape=tuple((d,n) for d,n,c in groups)
        free, images = determinant_packet(shape)
        projection[groups] = len(words),images
        for state in free:
            labels=[(d,c) for d,n,c in groups for _ in range(n)]
            words.append(tuple((d,c,a) for (d,c),a in zip(labels,state)))
    return words,projection

def boundary(words, projection, brackets):
    for word in words:
        values = Counter()
        for i,(d,c,a) in enumerate(word):
            for j in range(i+1,len(word)):
                e,b,t = word[j]
                if not brackets[c][b]:
                    continue
                rest=word[:i]+word[i+1:j]+word[j+1:]
                for g,n in brackets[c][b]:
                    merged=(d+e,g,a+t)
                    if merged in rest:
                        continue
                    groups,state,sign=packet((merged,)+rest)
                    offset,images=projection[groups]
                    mult=n*sign*(-1)**(i+j-1)
                    for k, coefficient in images[state].items():
                        values[offset+k] += mult*coefficient
        yield {k:n for k,n in values.items() if n}

def calculate(ce, weight, rank=rational_rank):
    l=ce.length
    low,left=basis(ce,l-1,weight)
    mid,middle=basis(ce,l,weight)
    high,_=basis(ce,l+1,weight)
    out=list(boundary(mid,left,ce.brackets))
    inc=list(boundary(high,middle,ce.brackets))
    for col in inc:
        check=Counter()
        for j,m in col.items():
            for i,n in out[j].items():
                check[i]+=m*n
        assert not any(check.values())
    ri,ro=rank(inc),rank(out)
    return dict(weight=weight,dimensions=[len(low),len(mid),len(high)],ranks=[ri,ro],homology=len(mid)-ri-ro)
