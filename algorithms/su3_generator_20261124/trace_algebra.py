"""Complete oriented relative-jet trace words for the B05 native grading.

No ansatz or reduced letter set is imposed: nonnegative flavour charges
(3,0,0) force every allowed letter to contain either no theta or theta_1.
The integer derivative-raising kernel imposes dotted spin zero.
"""
from collections import Counter, defaultdict
from functools import cache
from itertools import product
from math import gcd
from pathlib import Path
import sys

VENDOR = Path(__file__).resolve().parents[1] / "su3_independent_prime_20261118/runtime/work/cluster-next-v2.nosync/work/so8-bmn/vendor"
sys.path.insert(0, str(VENDOR))
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm
from fortuity.full.wordenum import multitraces_with_content
from flint import fmpz_mat


@cache
def contents(length, plus=4, minus=4, scalars=3):
    alphabet = sorted(((a,b,(0,) if s else ()) for a in range(plus+1)
                       for b in range(minus+1) for s in range(2)
                       if a+b+s and a+b+s <= plus+minus+scalars-length+1), reverse=True)
    output = []
    def visit(start, left, remaining, word):
        if not left:
            if remaining == (0,0,0): output.append(tuple(sorted(word)))
            return
        for index in range(start,len(alphabet)):
            letter=alphabet[index]; weight=letter[:2]+(len(letter[2]),)
            maximum=min([left]+[r//w for r,w in zip(remaining,weight) if w])
            for count in range(1,maximum+1):
                visit(index+1,left-count,tuple(r-count*w for r,w in zip(remaining,weight)),word+(letter,)*count)
    visit(0,length,(plus,minus,scalars),())
    return tuple(output)


@cache
def canonical_multi(multi):
    coefficient=1; traces=[]
    for trace in multi:
        sign,word=fw.canon_word(trace)
        if not sign: return 0,None
        coefficient*=sign; traces.append(word)
    sign,word=fw.canon_multi(traces)
    return (sign*coefficient,word) if sign else (0,None)


@cache
def conjugate_word(multi):
    sign=1; reverse=[]
    for trace in multi:
        degree=sum(map(is_ferm,trace))
        sign*=(-1)**(len(trace)+degree*(degree-1)//2)
        reverse.append(trace[::-1])
    cyclic,word=canonical_multi(tuple(reverse))
    return sign*cyclic,word


def conjugate(row):
    answer=Counter()
    for word,c in row.items():
        sign,other=conjugate_word(word)
        if sign: answer[other]+=sign*c
    return {w:c for w,c in answer.items() if c}


def parity_basis(words,parity=-1):
    """Solve the coefficient equations on every conjugation orbit.

    One coefficient is the free variable on a two-word orbit; its partner
    is fixed by the eigenvalue equation. No linear map is applied to seeds.
    """
    remaining=set(words)
    for word in sorted(words):
        if word not in remaining: continue
        sign,other=conjugate_word(word)
        assert other in remaining and conjugate_word(other)==(sign,word)
        remaining.discard(word); remaining.discard(other)
        if other != word: yield {word:1,other:sign*parity}
        elif sign == parity: yield {word:1}
    assert not remaining


def shape(multi):
    def trace_shape(trace):
        labels=tuple((a+b,len(theta)) for a,b,theta in trace)
        return min(v[i:]+v[:i] for v in (labels,labels[::-1]) for i in range(len(v)))
    return tuple(sorted(map(trace_shape,multi)))


@cache
def raised_word(word):
    result=Counter()
    for i,trace in enumerate(word):
        for j,(a,b,theta) in enumerate(trace):
            if not b: continue
            changed=trace[:j]+((a+1,b-1,theta),)+trace[j+1:]
            sign,target=canonical_multi(word[:i]+(changed,)+word[i+1:])
            if sign: result[target]+=sign*(a+1)
    return tuple((w,c) for w,c in result.items() if c)


def raise_row(row):
    result=Counter()
    for word,c in row.items():
        for changed,v in raised_word(word): result[changed]+=c*v
    return {w:c for w,c in result.items() if c}


def combine(recipe,rows):
    result=Counter()
    for i,c in recipe.items():
        for w,v in rows[i].items(): result[w]+=c*v
    return {w:c for w,c in result.items() if c}


def highest_weights(length,parity=-1,progress=None):
    blocks=defaultdict(set)
    for content in contents(length):
        for word in multitraces_with_content(content): blocks[shape(word)].add(word)
    if progress: progress(dict(stage="enumerated",length=length,words=sum(map(len,blocks.values())),shapes=len(blocks)))
    result=[]
    for index,key in enumerate(sorted(blocks)):
        rows=list(parity_basis(blocks[key],parity)); images=list(map(raise_row,rows))
        keys=sorted({w for row in images for w in row}); ids={w:i for i,w in enumerate(keys)}
        matrix=fmpz_mat(len(keys),len(rows))
        for j,row in enumerate(images):
            for w,c in row.items(): matrix[ids[w],j]=c
        kernel,dimension=matrix.nullspace()
        for j in range(dimension):
            recipe={i:int(kernel[i,j]) for i in range(len(rows)) if kernel[i,j]}
            divisor=gcd(*recipe.values()); recipe={i:c//divisor for i,c in recipe.items()}
            expr=combine(recipe,rows)
            assert not raise_row(expr)
            assert conjugate(expr)=={w:parity*c for w,c in expr.items()}
            result.append(expr)
        if progress and index%100==0: progress(dict(stage="highest_weights",length=length,blocks=index,total=len(blocks),singlets=len(result)))
    return result


def differential(row):
    return dict(fw.Q_expr(row))


def skeleton(row):
    words=iter(row); word=next(words)
    result=tuple(sorted(Counter((len(t),a+b) for trace in word for a,b,t in trace).items()))
    assert all(tuple(sorted(Counter((len(t),a+b) for trace in w for a,b,t in trace).items()))==result for w in words)
    return result
