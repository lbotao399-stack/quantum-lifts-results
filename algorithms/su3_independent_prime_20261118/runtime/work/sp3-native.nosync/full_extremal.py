"""Extremal projector on full letters in the existing Taylor convention.

Tolstoy, math-ph/0412087, (3.35)-(3.37), for SU(3), followed by the
commuting dotted SU(2) projector. This changes no Q coefficients.
"""
import sp_environment
import sys
sys.path.insert(0,str(sp_environment.WORK/'cluster-next-v2.nosync/work/so8-protected'))
from full_hw_so import so_canon
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache


@lru_cache(maxsize=16384)
def word_action(word,row,column):
    out=defaultdict(int)
    for trace_index,trace in enumerate(word):
        for letter_index,(a,b,theta) in enumerate(trace):
            if row=='spin':
                if column==1:
                    if not b:
                        continue
                    factor=a+1
                    letter=(a+1,b-1,theta)
                else:
                    assert column==-1
                    if not a:
                        continue
                    factor=b+1
                    letter=(a-1,b+1,theta)
            else:
                assert row!=column and row in range(3) and column in range(3)
                if column not in theta or row in theta:
                    continue
                factor=(-1)**sum(min(row,column)<i<max(row,column) for i in theta)
                letter=(a,b,tuple(sorted(row if i==column else i for i in theta)))
            changed=word[:trace_index]+(trace[:letter_index]+(letter,)+trace[letter_index+1:],)+word[trace_index+1:]
            sign,target=so_canon(changed)
            if sign:
                out[target]+=sign*factor
    return {word:c for word,c in out.items() if c}


def action(expr,row,column,prime=None):
    out=defaultdict(int)
    for word,c in expr.items():
        for target,factor in word_action(word,row,column).items():
            out[target]+=c*factor
    if prime is None:
        return {w:c for w,c in out.items() if c}
    return {w:int(c)%prime for w,c in out.items() if c%prime}


def root_factor(expr,row,column,weight,prime=None):
    shift=1 if row=='spin' else column-row
    lowering=('spin',-1) if row=='spin' else (column,row)
    raised=expr
    result=dict(expr)
    coefficient=1
    for power in range(1,1000):
        raised=action(raised,row,column,prime)
        if not raised:
            return result
        denominator=power*(weight+shift+power)
        assert denominator and (prime is None or denominator%prime)
        coefficient*=(-Fraction(1,denominator) if prime is None else -pow(denominator,-1,prime))
        if prime is not None:
            coefficient%=prime
        lowered=raised
        for _ in range(power):
            lowered=action(lowered,*lowering,prime)
        for word,c in lowered.items():
            result[word]=result.get(word,0)+coefficient*c
        result=({w:c for w,c in result.items() if c} if prime is None
                else {w:int(c)%prime for w,c in result.items() if c%prime})
    raise AssertionError('Nonterminating root string')


def project(expr,sector,prime=None):
    assert sector[0]>=sector[1] and sector[2]>=sector[3]>=sector[4]
    for row,column in ((1,2),(0,2),(0,1)):
        difference=sector[row+2]-sector[column+2]
        assert difference%2==0
        expr=root_factor(expr,row,column,difference//2,prime)
    difference=sector[0]-sector[1]
    assert difference%2==0
    return root_factor(expr,'spin',1,difference//2,prime)
