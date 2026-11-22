"""Invert the tree charge on orthogonal first-jet words containing an odd trace."""
from collections import Counter
from fractions import Fraction
from itertools import product

import mixed_incoming as traces
from fortuity.full.letters import Q_letter, is_ferm


def primitive(word):
    if not any(len(trace)%2 for trace in word):
        raise ValueError('An odd-length trace is required')
    if any(i+j+len(theta)!=1 for trace in word for i,j,theta in trace):
        raise ValueError('Every factor must be a first jet')
    index = next(i for i, trace in enumerate(word) if len(trace)%2)
    trace = list(word[index])
    previous = sum(is_ferm(letter) for t in word[:index] for letter in t)
    parity = 1
    result = Counter()
    for end in range(len(trace)-1,0,-1):
        for i in range(end):
            x,y=trace[i:i+2]
            swap = -1 if is_ferm(x) and is_ferm(y) else 1
            if x!=y or swap==-1:
                jet=(x[0]+y[0],x[1]+y[1],tuple(sorted(x[2]+y[2])))
                coefficients=Counter()
                for value,left,right in Q_letter(jet):
                    coefficients[left,right]+=value
                expected=Counter({(x,y):1})
                expected[y,x]-=swap
                factor=Fraction(expected[x,y],coefficients[x,y])
                assert all(factor*coefficients[k]==expected[k] for k in coefficients.keys()|expected.keys())
                sign=(-1)**(previous+sum(map(is_ferm,trace[:i])))
                replaced=tuple(trace[:i]+[jet]+trace[i+2:])
                orientation,term=traces.so_canon(word[:index]+(replaced,)+word[index+1:])
                if orientation:
                    result[term]+=Fraction(orientation*parity*sign,2)*factor
            trace[i],trace[i+1]=y,x
            parity*=swap
    return {w:n for w,n in result.items() if n}


if __name__ == '__main__':
    letters=((1,0,()),(0,1,()),(0,0,(0,)),(0,0,(1,)),(0,0,(2,)))
    def check(word):
        orientation,w=traces.so_canon(word)
        expected={w:orientation} if orientation else {}
        assert traces.q(primitive(word))==expected
        return int(bool(orientation))
    tested=sum(check((word,)) for word in product(letters,repeat=5))
    # Both source grades of interest, with independent spin and all three flavours.
    from random import Random
    shuffle=Random(512349).shuffle
    for power in (9,10):
        pool=[g for g in letters[:2] for _ in range(power-3)]+[g for g in letters[2:] for _ in range(3)]
        for _ in range(12):
            shuffle(pool)
            tested+=check((tuple(pool[:2]),tuple(pool[2:7]),tuple(pool[7:])))
            tested+=check((tuple(pool),))
    print(tested, 'nonzero words with exactly verified tree primitives')
