"""Evaluate Q1 on the collected F4 source at the archive's two families of colour components."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import cache, lru_cache
from hashlib import sha256
from itertools import product
import json
from multiprocessing import get_context
from pathlib import Path
from zipfile import ZipFile

from native import matrix, multiply, odd
from verify import ARCHIVE, contents

PSI = (0,0,-1)
PHI = (0,0,-2)


def marked_graphs(rows):
    """Move the two inputs to the front, then collect spectator traces before marked traces."""
    result = Counter()
    for row in rows:
        if row['profile'] != 'Psi':
            continue
        factors = [(kind,tuple(map(tuple,word))) for kind,word in row['factors']]
        leaves = [(i,j,jet) for i,(_,word) in enumerate(factors) for j,jet in enumerate(word)]
        psi_trace,psi_site,psi = next(x for x in leaves if x[2][2].bit_count() == 2)
        before = sum(odd(jet) for i,j,jet in leaves if (i,j)<(psi_trace,psi_site))
        # Positive one-loop unit, ordered psi_12 phi_3 and psi_23 phi_1 both have negative sign.
        orientation = 1 if psi[2] == 5 else -1
        for phi_trace,phi_site,phi in leaves:
            if phi != (0,0,7^psi[2]):
                continue
            normal, marked = [], []
            for i,(kind,word) in enumerate(factors):
                slots = list(word)
                if i == psi_trace:
                    slots[psi_site] = PSI
                if i == phi_trace:
                    slots[phi_site] = PHI
                parity = sum(odd(jet) for jet in slots if jet not in (PSI,PHI))
                group = marked if i in (psi_trace,phi_trace) else normal
                group.append((i,(kind,tuple(slots)),parity))
            exchanges = sum(p*q for i,_,p in marked for j,_,q in normal if j>i)
            key = tuple(x[1] for x in normal),tuple(x[1] for x in marked)
            result[key] += orientation*(-1)**(before+exchanges)*Fraction(row['coefficient'])
    return [(w,c) for w,c in result.items() if c]


class FirstAction:
    def __init__(self, tables, rows):
        self.dual = list(map(matrix,tables['nativeDualTranspose26']))
        self.metric = matrix(tables['inverseMetric'])
        self.ad = list(map(matrix,tables['adjointMatrices']))
        self.graphs = marked_graphs(rows)

    @lru_cache(maxsize=24000)
    def matrix_word(self, colours):
        if len(colours) == 1:
            return self.dual[colours[0]]
        return multiply(self.matrix_word(colours[:-1]),self.dual[colours[-1]])

    @lru_cache(maxsize=120000)
    def factor_value(self, kind, colours):
        if not kind:
            return self.metric.get(colours,0)
        return sum(v for (i,j),v in self.matrix_word(colours).items() if i == j)

    @cache
    def omega(self, u, v):
        row = Counter()
        for a,b in ((u,v),(v,u)):
            for (r,s),g in self.metric.items():
                for (i,h),x in self.ad[r].items():
                    if a == i:
                        for (j,k),y in self.ad[s].items():
                            if j == b:
                                row[h,k] += x*y*g/2
        return [(h,k,c) for (h,k),c in row.items() if c]

    def coefficient(self, indices):
        plus,minus,scalars = indices[:5],indices[5:10],indices[10:]
        assert all(tuple(sorted(set(v))) == tuple(v) for v in (plus,minus,scalars))
        leaves = (tuple(((1,0,0),u) for u in plus) + tuple(((0,1,0),v) for v in minus) +
                  tuple(((0,0,1<<i),t) for i,t in enumerate(scalars)))
        full = (1<<13)-1

        @cache
        def sign(a,b):
            ghosts = b & 1023
            return (-1)**sum((ghosts % (1<<i)).bit_count() for i in range(10) if a & (1<<i))

        def assignments(jet, allowed):
            return [(1<<i,c) for i,(s,c) in enumerate(leaves) if s == jet and (1<<i)&allowed]

        def ordered(choices):
            used = 0
            orientation = 1
            for bit,_ in choices:
                if bit & used:
                    return 0,0
                orientation *= sign(used,bit)
                used |= bit
            return used,orientation

        @cache
        def factor(f):
            kind,word = f
            values = Counter()
            for choices in product(*(assignments(j,full) for j in word)):
                used,s = ordered(choices)
                if s:
                    values[used] += s*self.factor_value(kind,tuple(c for _,c in choices))
            return {m:v for m,v in values.items() if v}

        @cache
        def spectators(normal):
            if not normal:
                return {0:1}
            values = Counter()
            for m,x in spectators(normal[:-1]).items():
                for n,y in factor(normal[-1]).items():
                    if not m & n:
                        values[m|n] += x*y*sign(m,n)
            return {m:x for m,x in values.items() if x}

        @cache
        def active_colour(traces,u,v):
            total = Fraction()
            for a,b,coefficient in self.omega(u,v):
                for kind,word in traces:
                    colours = tuple(a if k == -1 else b if k == -2 else k for k in word)
                    coefficient *= self.factor_value(kind,colours)
                    if not coefficient:
                        break
                total += coefficient
            return total

        @cache
        def active(marked,available,u,v):
            jets = [j for _,word in marked for j in word]
            options = [[(0,j[2])] if j in (PSI,PHI) else assignments(j,available) for j in jets]
            total = Fraction()
            for choices in product(*options):
                used,s = ordered(choices)
                if s:
                    assert used == available
                    colours = iter(c for _,c in choices)
                    traces = tuple((kind,tuple(next(colours) for _ in word)) for kind,word in marked)
                    total += s*active_colour(traces,u,v)
            return total

        total = Fraction()
        for (normal,marked),coefficient in self.graphs:
            for m,x in spectators(normal).items():
                for i,u in enumerate(plus):
                    if (1<<i)&m:
                        continue
                    for j,v in enumerate(minus,5):
                        emitted = (1<<i)|(1<<j)
                        if emitted & m or not self.omega(u,v):
                            continue
                        rest = full^(emitted|m)
                        total += coefficient*x*sign(emitted,m)*sign(emitted|m,rest)*active(marked,rest,u,v)
        return total


def inputs(archive):
    def read(tail):
        return json.loads(contents(archive,tail))
    tables = read('collected-O0/color-basis-native.json')
    rows = [json.loads(line) for line in contents(archive,'collected-O0/O0-canonical.jsonl').splitlines()]
    families = [read('cochains/inputs/exact-raw-anchors.json')['anchors'],
                read('exact-score/eta0-complete-exact-anchors.json')]
    return tables,rows,families


def direct(indices):
    return str(ACTION.coefficient(indices))


def replay(tables,rows,families,workers):
    global ACTION
    ACTION = FirstAction(tables,rows)
    # The archive stores six times a coefficient in its prescribed determinant orientation.
    expected = {}
    for family in families:
        for key,value in family:
            key = tuple(key)
            number = Fraction(value)
            assert expected.setdefault(key,number) == number
    keys = sorted(expected)
    with ProcessPoolExecutor(workers,mp_context=get_context('fork')) as pool:
        values = list(pool.map(direct,keys,chunksize=8))
    discrepancies = []
    for key,value in zip(keys,values):
        if 6*Fraction(value) != -expected[key]:
            discrepancies.append([key,value,str(expected[key])])
    return dict(family_sizes=list(map(len,families)),distinct_keys=len(keys),
                nonzero=sum(bool(Fraction(v)) for v in values),ratio_to_archived_anchor_over_six='-1',
                value_sha256=sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest(),
                discrepancies=discrepancies)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    args = parser.parse_args()
    if not __debug__ or args.workers<1:
        parser.error('positive worker count and assertions required')
    assert sha256(args.archive.read_bytes()).hexdigest() == ARCHIVE
    with ZipFile(args.archive) as z:
        data = inputs(z)
    print(json.dumps(replay(*data,args.workers),indent=2))
