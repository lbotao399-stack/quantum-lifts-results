"""Verify the lifted trace singlets before imposing any ghost commutation or finite trace identity."""
from collections import Counter
from functools import cache
import gzip
from itertools import product
import json
from math import prod

from flag_boundaries import HERE
from ghost_trace_quotient import sign
from trace_source_cycles import GLUE, PAIRS, canonical, odd
from trace_source_score import section


@cache
def raw_factor(w):
    par = sum(map(odd,w))
    possibilities = []
    for word,phase in ((w,1),(w[::-1],(-1)**(len(w)+par*(par-1)//2))):
        for i in range(len(w)):
            before = sum(map(odd,word[:i]))
            possibilities.append((word[i:]+word[:i],phase*(-1)**(before*(par-before))))
    target = min(k for k,_ in possibilities)
    signs = {s for k,s in possibilities if k==target}
    return (signs.pop(),target) if len(signs)==1 else (0,())


@cache
def raw_word(m):
    terms = [raw_factor(w) for w in m]
    factors = tuple(w for _,w in terms)
    odds = tuple(w for w in factors if sum(map(odd,w))%2)
    if len(set(odds))!=len(odds):
        return 0,()
    return prod(s for s,_ in terms)*sign(odds),tuple(sorted(factors))


def lift(m):
    out = Counter()
    for term in product(*(section(f) if f!=GLUE else [(GLUE,1)] for f in m)):
        s,key = raw_word(tuple(f for f,_ in term))
        if s:
            out[key] += s*prod(c for _,c in term)
    result = {k:x for k,x in out.items() if x}
    reduction = Counter()
    for k,x in result.items():
        s,w = canonical(k)
        if s:
            reduction[w] += s*x
    assert {k:v for k,v in reduction.items() if v}=={m:1}
    return result


@cache
def raise_word(m,a,b):
    out = Counter()
    for i,w in enumerate(m):
        for j,v in enumerate(w):
            if v==b:
                value,z = a,1
            elif a<3 and v>=5 and b in PAIRS[v-5] and a not in PAIRS[v-5]:
                labels = tuple(a if x==b else x for x in PAIRS[v-5])
                value,z = 5+PAIRS.index(tuple(sorted(labels))),sign(labels)
            else:
                continue
            s,key = raw_word(m[:i]+(w[:j]+(value,)+w[j+1:],)+m[i+1:])
            if s:
                out[key] += s*z
    return {k:x for k,x in out.items() if x}


def calculate():
    data = json.loads(gzip.decompress((HERE/"trace_source_cycles.json.gz").read_bytes()))
    sources = [lift(tuple(map(tuple,m))) for m in data["monomials"]]
    sizes = []
    for c in data["singlet_columns"]:
        raw = Counter()
        for i,x in c:
            for m,y in sources[i].items():
                raw[m] += x*y
        actual = {m:y for m,y in raw.items() if y}
        for a,b in ((3,4),(0,1),(1,2)):
            out = Counter()
            for m,y in actual.items():
                for k,z in raise_word(m,a,b).items():
                    out[k] += y*z
            assert not any(out.values())
        sizes.append(len(actual))
    return dict(monomials=len(sources),actual_singlet_supports=sizes,
                raw_raising_identities=3*len(sizes))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
