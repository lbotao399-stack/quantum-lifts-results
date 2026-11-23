"""Highest-weight trace polynomials for the F4 fifth-glueball partner calculation."""
from collections import defaultdict
from functools import cache
from itertools import combinations, permutations, product
from math import gcd
import numpy as np

HEIGHTS = (4, 4, 3)
MARK = 1 << 24
ADJOINT = 1 << 27
PAIRS = tuple(combinations(range(4), 2))


@cache
def flags():
    orders = [list(permutations(range(h))) for h in HEIGHTS]
    digits, signs = [], []
    for columns in product(*orders):
        digits.append(tuple(i for col in columns for i in col))
        inv = sum(a>b for col in columns for j,a in enumerate(col) for b in col[j+1:])
        signs.append((-1)**inv)
    return np.asarray(digits, dtype=np.int32).T.copy(), np.asarray(signs, dtype=np.int32)


def encode(ids):
    result = np.zeros(flags()[1].shape, dtype=np.int32)
    for i in ids:
        result = 4*result + flags()[0][i]
    return result


def rotmin(slots):
    n = len(slots)
    code = encode(slots)
    result = code.copy()
    for _ in range(1,n):
        code = (4*code) % 4**n + (code >> (2*n-2))
        result = np.minimum(result, code)
    return result


@cache
def factor(slots, adjoint):
    dots = [i for i,t in enumerate(slots) if len(t)==2]
    phase = np.ones(flags()[1].shape, dtype=np.int32)
    if dots:
        position, = dots
        i,j = slots[position]
        digits = flags()[0]
        phase *= np.sign(digits[j]-digits[i])
        lo = np.minimum(digits[i],digits[j])
        hi = np.maximum(digits[i],digits[j])
        tags = np.zeros(phase.shape, dtype=np.int32)
        for k,(a,b) in enumerate(PAIRS,1):
            tags += k*((lo==a)&(hi==b))
        raw = tuple(v[0] for v in slots[position+1:]+slots[:position])
        u,v = encode(raw), encode(raw[::-1])
        n = len(raw)
    else:
        raw = tuple(t[0] for t in slots)
        u,v = rotmin(raw),rotmin(raw[::-1])
        n = len(raw)
        tags = 0
    if len(slots)%2:
        phase *= np.sign(v-u)
    return (np.minimum(u,v)+4**n+MARK*tags+ADJOINT*adjoint, phase)


def partitions(n, limit=None):
    if n==0:
        yield ()
        return
    for j in range(2, min(limit or n,n)+1):
        for r in partitions(n-j,j):
            yield (j,)+r


def raw_specs(mode, seed, two_reps):
    rng=np.random.default_rng(731921+seed)
    if mode=='out':
        kind='scalar'
        y=11
    else:
        kind=('two','two','three')[seed%3]
        y=9
    shapes=tuple(partitions(y))
    shape=shapes[(seed if mode=='out' else seed//3)%len(shapes)]
    order=tuple(int(i) for i in rng.permutation(11))
    if kind=='scalar':terms=[(1,[(x,) for x in order])]
    elif kind=='three':
        a,b,c=order[:3]
        tokens=[(a,b,c)]+[(i,) for i in order[3:]]
    else:
        tokens=[tuple(order[:2]),tuple(order[2:4])]+[(i,) for i in order[4:]]
    if kind!='scalar':
        rng.shuffle(tokens)
        if kind=='three':
            i=next(k for k,t in enumerate(tokens) if len(t)==3)
            a,b,c=tokens[i]
            terms=[]
            for x,y0,z in ((a,b,c),(b,c,a),(c,a,b)):
                terms += [(1, tokens[:i]+[(x,),(y0,z)]+tokens[i+1:]),
                          (-1,tokens[:i]+[(y0,z),(x,)]+tokens[i+1:])]
            expanded_positions=[i]*len(terms)
        else:
            terms=[]; expanded_positions=[]
            for number,i in enumerate(k for k,t in enumerate(tokens) if len(t)==2):
                a,b=tokens[i]
                for s,t in ((1,[(a,),(b,)]),(-1,[(b,),(a,)])):
                    terms.append((s*(-1)**number,tokens[:i]+t+tokens[i+1:]))
                    expanded_positions.append(i)
    reps=rng.integers(0,2,len(shape)) if two_reps else [0]*len(shape)
    for term,(sign,letters) in enumerate(terms):
        factors=[];past=0;begin=0
        for sz,rep in zip(shape,reps):
            n=sz+(kind!='scalar' and past<=expanded_positions[term]<past+sz)
            factors.append((tuple(letters[begin:begin+n]),int(rep)))
            begin+=n;past+=sz
        assert begin==len(letters)
        yield sign,factors


def polynomial(mode,seed,two_reps):
    keys=[];coeff=[]
    for sign,specs in raw_specs(mode,seed,two_reps):
        rows=[]
        phase=flags()[1].copy()*sign
        for slots,rep in specs:
            ids,sg=factor(slots,rep)
            rows.append(ids);phase*=sg
        active=phase!=0
        if np.any(active):
            keys.append(np.sort(np.stack(rows,axis=1)[active],axis=1))
            coeff.append(phase[active])
    factor.cache_clear()
    if not keys:return {}
    order,ix=np.unique(np.concatenate(keys),axis=0,return_inverse=True)
    co=np.zeros(len(order),dtype=np.int64)
    np.add.at(co,ix,np.concatenate(coeff))
    good=np.flatnonzero(co)
    if not len(good):return {}
    content=gcd(*(int(co[i]) for i in good))
    return {tuple(map(int,order[i])):int(co[i])//content for i in good}
