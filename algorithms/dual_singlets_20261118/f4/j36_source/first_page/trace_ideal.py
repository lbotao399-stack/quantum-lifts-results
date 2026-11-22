"""Exact graded trace quotient after imposing the polarized F4 quartic on Lie words."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import permutations, product
import json
from math import gcd

from ghost_trace_quotient import rotate, sign, trace_counts

ZERO = (0, 0, 0, 0, 0)


@cache
def trace(word, commute=False):
    """Move the two odd matrices to the front by the ghost-splitting relations."""
    ghosts = tuple(k for k in word if k >= 3)
    if len(ghosts) != len(set(ghosts)) or len(word) < 2:
        return 0, ()
    scalars = tuple(k for k in word if k < 3)
    if commute:
        scalars = tuple(sorted(scalars))
        forward = backward = scalars
    else:
        forward, backward = rotate(scalars), rotate(scalars[::-1])
    reverse_sign = (-1)**len(word)
    if forward == backward and reverse_sign < 0:
        return 0, ()
    phase = sign(ghosts) * (1 if forward <= backward else reverse_sign)
    return phase, tuple(sorted(ghosts)) + min(forward, backward)


@cache
def parity(factor):
    return sum(k >= 3 for k in factor) % 2


@cache
def order(factors):
    odd = tuple(f for f in factors if parity(f))
    if len(odd) != len(set(odd)):
        return 0, ()
    return sign(odd), tuple(sorted(factors))


def degree(word):
    return tuple(word.count(i) for i in (3, 4, 0, 1, 2))


def difference(left, right):
    return tuple(a-b for a, b in zip(left, right))


@cache
def generators():
    factors = set()
    for k in range(7):
        for w in product(range(3), repeat=k):
            for ghosts in ((), (3,), (4,), (3,4)):
                s, f = trace(ghosts+w)
                if s:
                    factors.add(f)
    return tuple((f, degree(f)) for f in sorted(factors))


@cache
def monomials(weight):
    if weight == ZERO:
        return ((),)
    out = []
    for f, d in generators():
        rest = difference(weight, d)
        if min(rest) < 0:
            continue
        for tail in monomials(rest):
            if tail and (tail[0] < f or (tail[0] == f and parity(f))):
                continue
            out.append((f,)+tail)
    return tuple(out)


def lyndon(w):
    return all(w < w[i:] for i in range(1,len(w)))


@cache
def lie(w):
    if len(w) == 1:
        return {w:1}
    i = next(j for j in range(1,len(w)) if lyndon(w[j:]))
    out = Counter()
    for u, a in lie(w[:i]).items():
        for v, b in lie(w[i:]).items():
            out[u+v] += a*b
            out[v+u] -= a*b
    return {k:v for k,v in out.items() if v}


def trace_product(arguments, groups):
    ids = tuple(i for group in groups for i in group)
    odd = tuple(i for i in ids if arguments[i][0] >= 3)
    phase = sign(odd)
    expansion = {():phase}
    for group in groups:
        terms = Counter()
        for factors in product(*(lie(arguments[i]).items() for i in group)):
            w = ()
            coefficient = 1
            for word, v in factors:
                w += word
                coefficient *= v
            s, f = trace(w)
            if s:
                terms[f] += s*coefficient
        enlarged = Counter()
        for a, x in expansion.items():
            for f, y in terms.items():
                s, key = order((*a,f))
                if s and y:
                    enlarged[key] += s*x*y
        expansion = enlarged
    return expansion


def normalized(row):
    items = sorted((k,v) for k,v in row.items() if v)
    if not items:
        return ()
    divisor = gcd(*(v for _,v in items)) * (1 if items[0][1] > 0 else -1)
    return tuple((k,v//divisor) for k,v in items)


@cache
def relations(ratio):
    # The longest factor uses five scalar letters when the other three
    # arguments are a scalar and the two distinct dotted ghosts.
    words = [w for n in range(1,6) for w in product(range(3),repeat=n) if lyndon(w)]
    def choose(n, room, start):
        if not n:
            yield ()
            return
        for j in range(start,len(words)):
            w = words[j]
            if n*len(w) > room:
                break
            for tail in choose(n-1,room-len(w),j):
                yield (w,)+tail
    by_weight = defaultdict(set)
    for ghosts in ((),((3,),),((4,),),((3,),(4,))):
        for args in choose(4-len(ghosts),6,0):
            args += ghosts
            wt = degree(sum(args, ()))
            if max(wt[2:]) > 4:
                continue
            relation = Counter()
            for p in permutations((1,2,3)):
                for key,v in trace_product(args, ((0,*p),)).items():
                    relation[key] += ratio.denominator*v
            for a,b,c in ((1,2,3),(2,1,3),(3,1,2)):
                for key,v in trace_product(args, ((0,a),(b,c))).items():
                    relation[key] -= 2*ratio.numerator*v
            packed = normalized(relation)
            if packed:
                by_weight[wt].add(packed)
    return by_weight


def echelon(basis, row):
    row = {k:v for k,v in row.items() if v}
    while row:
        j = min(row)
        x = row[j]
        if j not in basis:
            basis[j] = {k:Fraction(v,x) for k,v in row.items()}
            return
        for k,v in basis[j].items():
            y = row.get(k,0)-x*v
            if y:
                row[k] = y
            else:
                row.pop(k,None)


@cache
def abelian(monomial):
    phase = 1
    factors = []
    for f in monomial:
        s, factor = trace(f,True)
        if not s:
            return 0, ()
        phase *= s
        factors.append(factor)
    s, m = order(tuple(factors))
    return s*phase, m


def calculation(weight, ratio):
    source = monomials(weight)
    index = {m:i for i,m in enumerate(source)}
    assert len(index) == len(source) == sum(trace_counts(False)[(*weight,k)] for k in range(1,8))
    ideal = {}
    for wt, rows in relations(ratio).items():
        rest = difference(weight,wt)
        if min(rest) < 0:
            continue
        for spectator in monomials(rest):
            for relation in rows:
                row = Counter()
                for m, v in relation:
                    s, key = order(m+spectator)
                    if s:
                        row[index[key]] += s*v
                echelon(ideal,row)
    initial = len(ideal)
    section = {}
    by_traces = {}
    for k in range(1,8):
        for j,m in enumerate(source):
            if len(m) != k:
                continue
            s, key = abelian(m)
            if not s:
                echelon(ideal,{j:1})
            elif key in section:
                previous, orientation = section[key]
                echelon(ideal,{j:1,previous:-s*orientation})
            else:
                section[key] = (j,s)
        by_traces[k] = len(ideal)-initial
    assert len(section) == sum(trace_counts(True)[(*weight,k)] for k in range(1,8))
    return dict(weight=weight, formal=len(source), quartic_rank=initial,
                scalar_commutator_images_up_to_each_trace_count=list(by_traces.values()))


def calculate():
    outputs = []
    for dim,ratio in ((26,Fraction(1,12)),(52,Fraction(5,108))):
        rows, totals = [], [0]*7
        for shift in (0,1):
            for perm in permutations(range(3)):
                weight = (4+shift,4-shift,*(2+v-i for i,v in enumerate(perm)))
                row = calculation(weight,ratio)
                rows.append(row)
                for i,v in enumerate(row["scalar_commutator_images_up_to_each_trace_count"]):
                    totals[i] += (-1)**shift*sign(perm)*v
        assert totals == sorted(totals)
        outputs.append(dict(representation=dim,trace_ratio=str(ratio),weights=rows,
                            singlet_images_up_to_each_trace_count=totals))
    return outputs


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
