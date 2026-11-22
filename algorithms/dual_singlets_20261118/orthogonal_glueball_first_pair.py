"""A first quantum pairing on the extremal orthogonal determinant source."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import factorial, prod
from pathlib import Path

import mixed_incoming as traces
from orthogonal_odd_trace_primitive import primitive
from rank_four_lowered_replay import Q1d_multi, terminal
from terminal_commutators import sign


PHI = tuple((0,0,(i,)) for i in range(3))
PLUS, MINUS = (1,0,()), (0,1,())
GHOSTS = (PLUS, MINUS)
# H, E, F in the fundamental representation of the auxiliary SL(2).
ENTRIES = (((0,0,1),(1,1,-1)), ((0,1,1),), ((1,0,1),))


def scalar_code(indices):
    return sum(1 << (3*i) for i in indices)


class PfaffianCoefficient:
    """Coefficients of Pf(J+A)/Pf(J) in ordered colour traces."""

    def __init__(self, order, doubles):
        self.doubles = doubles
        det = [(tuple(3*i+j for i,j in enumerate(p)), sign(p)) for p in permutations(range(3))]
        powers = Counter()
        for terms in product(det, repeat=order):
            indices = tuple(i for block, _ in terms for i in block)
            powers[tuple(sorted(indices))] += prod(s for _,s in terms)
        self.weights = {}
        self.divisors = set()
        for indices, coefficient in powers.items():
            if not coefficient:
                continue
            self.weights[scalar_code(indices)] = coefficient*prod(factorial(n) for n in Counter(indices).values())
            partial = {0}
            for i in indices:
                partial.update([x+(1 << (3*i)) for x in partial])
            self.divisors.update(partial)

    @lru_cache(None)
    def exterior_sign(self, left, right):
        return (-1)**sum((left >> (j+1)).bit_count() for j in range(6) if right & (1 << j))

    def multiply(self, left, right):
        result = Counter()
        for (m,x,z), a in left.items():
            for (n,y,v), b in right.items():
                if m & n or z+v > self.doubles or x+y not in self.divisors:
                    continue
                result[m|n,x+y,z+v] += a*b*self.exterior_sign(m,n)
        return {k:v for k,v in result.items() if v}

    @lru_cache(None)
    def aux_trace(self, tiles):
        # Each nonnegative entry is a single colour letter. A double ghost
        # path has identity auxiliary matrix, with its orientation recorded
        # separately by tilings().
        states = {(0,0,0,0):1, (1,1,0,0):1}
        for letter in tiles:
            new = Counter()
            for (start,end,mask,weight), value in states.items():
                for s, entries in enumerate(ENTRIES):
                    bit = 1 << (3*(letter-3)+s) if letter >= 3 else 0
                    if mask & bit:
                        continue
                    degree = weight + ((1 << (3*(3*letter+s))) if letter < 3 else 0)
                    if degree not in self.divisors:
                        continue
                    for i,j,c in entries:
                        if end==i:
                            new[start,j,mask|bit,degree] += c*value*self.exterior_sign(mask,bit)
            states = {k:v for k,v in new.items() if v}
        result = Counter()
        for (i,j,m,x), value in states.items():
            if i==j:
                result[m,x] += value
        return {k:v for k,v in result.items() if v}

    @lru_cache(None)
    def tilings(self, letters):
        if not letters:
            return [((),0,1)]
        rows = [((letters[0],)+seq,z,c) for seq,z,c in self.tilings(letters[1:])]
        if letters[:2] in ((3,4),(4,3)):
            orientation = 1 if letters[0]==3 else -1
            rows.extend((seq,z+1,c*orientation) for seq,z,c in self.tilings(letters[2:])
                        if z<self.doubles)
        return rows

    @lru_cache(None)
    def trace(self, word):
        # Expand 1/2 tr log(1+A), including every rotation and reversal
        # exactly once as an ordered word before canonicalizing its colour.
        orbit = {v[i:]+v[:i] for v in (word,word[::-1]) for i in range(len(word))}
        poly = Counter()
        for raw in orbit:
            orientation, canonical = traces.so_canon((raw,))
            if not orientation:
                continue
            assert canonical==(word,)
            encoded = tuple((PHI+GHOSTS).index(letter) for letter in raw)
            for seq,z,c in self.tilings(encoded):
                size = len(seq)+z
                factor = Fraction(orientation*c*(-1)**(size+1), 2*size)
                for (m,x),v in self.aux_trace(seq).items():
                    poly[m,x,z] += factor*v
        return {k:v for k,v in poly.items() if v}

    @lru_cache(None)
    def traces_product(self, words):
        if not words:
            return {(0,0,0):1}
        return self.multiply(self.traces_product(words[:-1]),self.trace(words[-1]))

    @lru_cache(None)
    def coefficient(self, words):
        # Each of the six auxiliary odd markers is written before its odd
        # colour matrix. Moving all markers to the left gives (-1)^15.
        multiplicity = prod(factorial(n) for n in Counter(words).values())
        poly = self.traces_product(words)
        value = sum(c*poly.get((63,x,self.doubles),0) for x,c in self.weights.items())
        return -Fraction(value,multiplicity)


def orientations(letters):
    out = set()
    for w in set(permutations(letters)):
        s, raw = traces.so_canon((w,))
        if s:
            out.add(raw[0])
    return sorted(out)


def remove_once(word, subset):
    result = list(word)
    for item in subset:
        result.remove(item)
    return tuple(result)


@lru_cache(None)
def single_quantum(colours, active):
    result = Counter()
    for w, c in primitive((active,)).items():
        for v, a in Q1d_multi(w,colours,'SO_odd').items():
            s, monomial = terminal(v)
            if s:
                result[monomial] += c*a*s
    return {k:v for k,v in result.items() if v}


def highest_product(colours, power, order):
    # For orders one and three the diagnostic has only even trace factors:
    # all are bilinears of ghosts, or pairs of equal scalars.
    target = ((PLUS,MINUS),)*(power-order+1)
    target += tuple((p,p) for p in PHI for _ in range((order-1)//2))
    sign_target, target = traces.so_canon(target)
    assert sign_target
    pf = PfaffianCoefficient(order,power-order-3)
    rows = []
    for size in (1,2):
        for subset in sorted(set(combinations(target,size))):
            leaves = tuple(x for w in subset for x in w)
            if PLUS not in leaves or MINUS not in leaves:
                continue
            spectator = remove_once(target, subset)
            # Spectators have even parity and are already canonically ordered.
            s, wanted = terminal(subset)
            for active in orientations(remove_once(leaves,GHOSTS)+PHI):
                q = s*single_quantum(colours,active).get(wanted,0)
                if not q:
                    continue
                source = tuple(sorted((active,)+spectator))
                coefficient = pf.coefficient(source)
                if coefficient:
                    rows.append(dict(outputs=size,active=[(PHI+GHOSTS).index(x) for x in active],
                                     relation=str(coefficient),quantum=str(q),
                                     contribution=str(q*coefficient)))
    return dict(colours=colours,power=power,order=order,
                target=[[(PHI+GHOSTS).index(x) for x in w] for w in target],
                input_letters=['phi1','phi2','phi3','plus','minus'],
                coefficient=str(sum(Fraction(r['contribution']) for r in rows)), terms=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--colours',type=int,required=True)
    p.add_argument('--power',type=int,required=True)
    p.add_argument('--order',type=int,choices=(1,3),required=True)
    p.add_argument('--out',type=Path,required=True)
    a = p.parse_args()
    if not __debug__ or a.colours != a.power+2*a.order+2 or a.power < a.order+3:
        p.error('Use the extremal orthogonal family with assertions enabled')
    result=highest_product(a.colours,a.power,a.order)
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(a.colours,a.power,a.order,len(result['terms']),result['coefficient'])
