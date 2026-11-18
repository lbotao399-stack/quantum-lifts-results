"""Graded cyclic SU trace algebra, with charge conjugation kept as an action.

The pinned SO/Sp modules are never modified.  An isolated instance of the
extremal projector is bound to cyclic canonicalization before any call.
"""
import sp_environment
from collections import defaultdict
from functools import lru_cache
import importlib.util
from pathlib import Path
from fortuity.full import words as fw
from fortuity.full.letters import is_ferm


@lru_cache(maxsize=32768)
def canonical_multi(multi):
    coefficient, traces = 1, []
    for trace in multi:
        sign,word = fw.canon_word(trace)
        if not sign:
            return 0,None
        coefficient *= sign
        traces.append(word)
    sign,word = fw.canon_multi(traces)
    return (sign*coefficient,word) if sign else (0,None)


def canonical(expr):
    result=defaultdict(int)
    for multi,coefficient in expr.items():
        sign,word=canonical_multi(multi)
        if sign:
            result[word] += sign*coefficient
    return {word:coefficient for word,coefficient in result.items() if coefficient}


@lru_cache(maxsize=32768)
def conjugate_word(multi):
    sign=1
    reverse=[]
    for trace in multi:
        degree=sum(is_ferm(letter) for letter in trace)
        sign *= (-1)**(len(trace)+degree*(degree-1)//2)
        reverse.append(trace[::-1])
    cyclic,word=canonical_multi(tuple(reverse))
    return sign*cyclic,word


def conjugate(expr):
    result=defaultdict(int)
    for word,coefficient in expr.items():
        sign,target=conjugate_word(word)
        if sign:
            result[target] += sign*coefficient
    return {word:coefficient for word,coefficient in result.items() if coefficient}


def parity_project(expr,parity):
    """Unnormalized (1 + parity*C); normalization of Q is unchanged."""
    assert parity in (-1,1)
    result=defaultdict(int,canonical(expr))
    for word,coefficient in conjugate(expr).items():
        result[word] += parity*coefficient
    return {word:coefficient for word,coefficient in result.items() if coefficient}


spec=importlib.util.spec_from_file_location('_su_cyclic_extremal_v1',Path(__file__).with_name('full_extremal.py'))
projection=importlib.util.module_from_spec(spec)
spec.loader.exec_module(projection)
projection.so_canon=canonical_multi


def differential(expr):
    return canonical(fw.Q_expr(expr))
