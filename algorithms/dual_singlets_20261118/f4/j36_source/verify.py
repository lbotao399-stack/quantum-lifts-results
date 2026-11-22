"""Differentiate the stored rational preimage and prove the required F4 trace identity."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import combinations, permutations, product
import gzip
import json
from pathlib import Path
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ARCHIVE = '0a06b306a12abbe77cadbcfab616c08c579d82086e9798aece27b8c83a0f7b52'


def parity(jet):
    return 1 - jet[2].bit_count() % 2


def koszul(order, odd):
    inversions = sum(i > j for k, i in enumerate(order) if odd[i]
                     for j in order[k+1:] if odd[j])
    return (-1)**inversions


@cache
def orient(trace):
    """Cyclicity and reversal for odd/even coefficients times orthogonal generators."""
    n = len(trace)
    odd = tuple(map(parity, trace))
    candidates = {}
    for direction in (1, -1):
        for start in range(n):
            order = tuple((start + direction*i) % n for i in range(n))
            key = tuple(trace[i] for i in order)
            sign = koszul(order, odd) * (1 if direction == 1 else (-1)**n)
            if key in candidates and sign != candidates[key]:
                return (), 0
            candidates[key] = sign
    key = min(candidates)
    return key, candidates[key]


def canonical(traces):
    keys = []
    odd = []
    sign = 1
    for trace in traces:
        t, s = orient(trace)
        if not s:
            return (), 0
        keys.append(t)
        sign *= s
        odd.append(sum(map(parity, t)) % 2)
    order = sorted(range(len(keys)), key=keys.__getitem__)
    if any(keys[i] == keys[j] and odd[i] for i,j in zip(order, order[1:])):
        return (), 0
    return tuple(keys[i] for i in order), sign*koszul(order, odd)


def insert(row, traces, coefficient):
    word, s = canonical(traces)
    if s:
        row[word] += s * coefficient


def cleaned(row):
    return {k: v for k, v in row.items() if v}


@cache
def replace_five(trace):
    if len(trace) in (2, 3):
        return { (trace,): Fraction(1) }
    assert len(trace) == 5
    row = Counter({(trace,): Fraction(4,3)})
    odd = tuple(map(parity, trace))
    for selection in range(32):
        left = [i for i in range(5) if selection & (1<<i)]
        right = [i for i in range(5) if not selection & (1<<i)]
        u, v = tuple(trace[i] for i in left), tuple(trace[i] for i in right)
        s = koszul(left + right, odd)
        row[(u+v,)] -= Fraction(s,12)
        if len(u) >= 2 and len(v) >= 2:
            row[(u,v)] += Fraction(s,72)
    return cleaned(row)


def source(rows):
    result = Counter()
    for row in rows:
        parts = []
        c = Fraction(row['coefficient'])
        for kind, raw in row['factors']:
            trace = tuple(map(tuple, raw))
            if kind == 0:
                assert len(trace) == 2
                c /= 6
            else:
                assert kind == 1
            parts.append(replace_five(trace).items())
        for selected in product(*parts):
            traces = []
            coefficient = c
            for words, v in selected:
                traces.extend(words)
                coefficient *= v
            insert(result, traces, coefficient)
    return cleaned(result)


def quartic_commutator(pair, labels):
    """Coefficient of four distinct polarizing parameters in 12 tr(X^4) - tr(X^2)^2."""
    i, j = pair
    singles = [(e,) for e in labels if e not in pair]
    row = Counter()
    for bracket, s in [((i,j), 1), ((j,i), -1)]:
        args = [bracket] + singles
        for a,b,c,d in permutations(args):
            insert(row, (a+b+c+d,), 12*s)
            insert(row, (a+b, c+d), -s)
    return cleaned(row)


def check_five():
    labels = tuple((i,0,1) for i in range(5))  # five even, independently labelled matrices
    difference = Counter()
    for terms, c in replace_five(labels).items():
        insert(difference, terms, c)
    insert(difference, (labels,), -1)
    pairs = list(combinations(labels, 2))
    rows = [quartic_commutator(pair, labels) for pair in pairs]
    columns = sorted(set(difference).union(*(r.keys() for r in rows)))
    system = [[Fraction(r.get(k,0)) for r in rows]+[difference[k]] for k in columns]
    pivots = []
    for j in range(len(rows)):
        i = next((i for i in range(len(pivots),len(system)) if system[i][j]), None)
        if i is None:
            continue
        pivot = len(pivots)
        system[i],system[pivot] = system[pivot],system[i]
        v = system[pivot][j]
        system[pivot] = [x/v for x in system[pivot]]
        for i,row in enumerate(system):
            if i != pivot and row[j]:
                v = row[j]
                system[i] = [x-v*y for x,y in zip(row,system[pivot])]
        pivots.append(j)
    coeff = [Fraction() for _ in rows]
    for i,j in enumerate(pivots):
        coeff[j] = system[i][-1]
    assert all(difference[k] == sum(c*r.get(k,0) for c,r in zip(coeff,rows)) for k in columns)
    return dict(identities=len(rows), rank=len(pivots), coordinates=len(columns),
                multipliers=list(map(str,coeff)))


@cache
def native_splits(jet):
    p, q, mask = jet
    result = []
    for left in range(8):
        if left & mask != left:
            continue
        right = mask ^ left
        inv = sum(i>j for i in range(3) if left & (1<<i)
                  for j in range(3) if right & (1<<j))
        e = mask.bit_count() + (1+left.bit_count())*right.bit_count() + inv
        for a,b in product(range(p+1), range(q+1)):
            u, v = (a,b,left), (p-a,q-b,right)
            if any(u) and any(v):
                # The minus sign is the anti-representation proved on the 52 generators.
                result.append((-((-1)**e), u, v))
    return result


def differentiate(rows):
    result = Counter()
    for c, traces in rows:
        before = 1
        for t,trace in enumerate(traces):
            for i,jet in enumerate(trace):
                for s, u, v in native_splits(jet):
                    new = trace[:i]+(u,v)+trace[i+1:]
                    insert(result, traces[:t]+(new,)+traces[t+1:], c*before*s)
                if parity(jet):
                    before *= -1
    return cleaned(result)


def contents(archive, tail):
    names = [f for f in archive.namelist() if f.endswith('/'+tail)]
    assert len(names) == 1
    return archive.read(names[0])


def load_primitive():
    return [(Fraction(c),tuple(tuple(map(tuple,t)) for t in w))
            for c,w in json.loads(gzip.decompress((HERE/'tree_preimage.json.gz').read_bytes()))]


def verify(archive):
    from native import check_colour, components
    assert sha256(archive.read_bytes()).hexdigest() == ARCHIVE
    primitive = load_primitive()
    with ZipFile(archive) as z:
        def read(name):
            return json.loads(contents(z,'collected-O0/'+name))
        raw = [json.loads(line) for line in contents(z,'collected-O0/O0-canonical.jsonl').splitlines()]
        matrices = read('color-basis-native.json')
        colour, dual, bracket = check_colour(matrices)
        pins = {name: read('example-'+name+'.json')['jets'] for name in ('Psi-2232','Dphi','D2c')}
    quintic = check_five()
    original = source(raw)
    image = differentiate(primitive)
    assert image == original
    assert not differentiate((c,w) for w,c in original.items())
    max_three = lambda word: sum(sorted((j[2].bit_count() for t in word for j in t),reverse=True)[:3])
    budget = [max(max_three(w) for c,w in primitive),max(max_three(w) for w in original)]
    assert budget == [5,4]
    assert all(sum(map(len,w)) == 12 and
               tuple(sum(j[t] for tr in w for j in tr) for t in (0,1)) == (4,4) and
               tuple(sum(bool(j[2] & (1<<i)) for tr in w for j in tr) for i in range(3)) == (2,2,2)
               for c,w in primitive if c)
    native = components(primitive, dual, bracket, pins)
    assert native == {'Psi-2232':'2232','Dphi':'-24','D2c':'5382'}
    return dict(colour=colour, five_letter_identity=quintic, supplied_terms=len(raw),
                tree_preimage_terms=len(primitive), canonical_source_words=len(original),
                residual_words=0, maximum_flavour_slots_on_three_inputs_X_O0=budget, native_components=native)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error('assertions must be enabled')
    print(json.dumps(verify(args.archive),indent=2))
