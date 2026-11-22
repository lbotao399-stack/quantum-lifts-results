"""Change from products of defining traces to monomials on distinct Cartan planes."""
from collections import Counter
from functools import lru_cache
from math import factorial, prod

from flint import fmpq_mat
from terminal_commutators import parity, words


@lru_cache(None)
def groupings(n):
    if not n:
        return ((),)
    result=[]
    for groups in groupings(n-1):
        result.append(groups+((n-1,),))
        for i,g in enumerate(groups):
            result.append(groups[:i]+(g+(n-1,),)+groups[i+1:])
    return tuple(result)


def matrix(power):
    weight=(power-2,power-2,2,2,2)
    basis=words(weight)
    indices={w:i for i,w in enumerate(basis)}
    rows=[]
    for w in basis:
        row=[0]*len(basis)
        factor=2**len(w)*(-1)**(sum(weight)//2)
        for partition in groupings(len(w)):
            blocks=[tuple(sum(w[i][j] for i in group) for j in range(5)) for group in partition]
            if any(max(block[:2])>1 for block in blocks):
                continue
            multiplicities=Counter(blocks)
            if any(n>1 and parity(block) for block,n in multiplicities.items()):
                continue
            order=sorted(range(len(blocks)),key=lambda i:(sum(blocks[i]),blocks[i]))
            cells={k:j for j,g in enumerate(order) for k in partition[g]}
            labels=[2*cells[i]+k for i,trace in enumerate(w) for k in range(2) if trace[k]]
            sign=(-1)**sum(i>j for n,i in enumerate(labels) for j in labels[n+1:])
            target=tuple(blocks[i] for i in order)
            row[indices[target]]+=factor*sign*prod(factorial(n) for n in multiplicities.values())
        assert row[indices[w]] and all(not v or len(basis[i])<=len(w) for i,v in enumerate(row))
        rows.append(row)
    result=fmpq_mat(rows)
    assert result.rank()==len(basis)
    return result
