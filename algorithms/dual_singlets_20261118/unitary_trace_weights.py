"""SU(3) highest weights of oriented trace polynomials at fixed charge conjugation."""
import math
from collections import Counter, defaultdict

from flint import fmpz_mat
import mixed_incoming as m
import su_trace_adapter_v1 as algebra


def parity_rows(words, parity):
    remaining=set(words)
    for word in sorted(words):
        if word not in remaining:
            continue
        sign,partner=algebra.conjugate_word(word)
        assert partner in remaining and algebra.conjugate_word(partner)==(sign,word)
        remaining.discard(word)
        remaining.discard(partner)
        if word!=partner:
            yield {word:1,partner:sign*parity}
        elif sign==parity:
            yield {word:1}
    assert not remaining


def highest_weights(weight, length, parity):
    blocks=defaultdict(set)
    for content in m.contents(weight,length):
        for word in m.multitraces_with_content(content):
            sign,canonical=algebra.canonical_multi(word)
            if sign:
                blocks[m.shape(canonical,"spin")].add(canonical)
    result=defaultdict(list)
    for block in blocks.values():
        sources=list(parity_rows(block,parity))
        if not sources:
            continue
        equations=[]
        for row in sources:
            terms={}
            for root in ((0,1),(1,2),("spin",1)):
                raised=algebra.projection.action(row,*root)
                terms.update({(root,w):n for w,n in raised.items()})
            equations.append(terms)
        keys={key for eq in equations for key in eq}
        ids={key:i for i,key in enumerate(sorted(keys,key=repr))}
        mat=fmpz_mat(len(ids),len(equations))
        for j,eq in enumerate(equations):
            for key,n in eq.items():
                mat[ids[key],j]=n
        ker,ncols=mat.nullspace()
        for j in range(ncols):
            coeff=[int(ker[i,j]) for i in range(len(sources))]
            gcd=math.gcd(*coeff)
            assert gcd and (mat*fmpz_mat([[v] for v in coeff])).is_zero()
            expr=Counter()
            for row,c in zip(sources,coeff):
                if c:
                    for w,n in row.items():
                        expr[w]+=n*(c//gcd)
            expr={w:n for w,n in expr.items() if n}
            assert algebra.conjugate(expr)=={w:parity*n for w,n in expr.items()}
            degrees={m.nferm(w) for w in expr}
            assert len(degrees)==1
            result[degrees.pop()].append(expr)
    return result
