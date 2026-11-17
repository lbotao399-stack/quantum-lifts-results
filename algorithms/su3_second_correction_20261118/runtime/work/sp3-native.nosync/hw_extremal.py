"""SU(3) extremal projection in the existing graded trace/Pfaffian algebra.

Formula: V. N. Tolstoy, 'Fortieth Anniversary of Extremal Projector Method',
eqs. (3.35)-(3.37), P=P12 P13 P23. Each sum terminates on these finite
polynomial representations. No change to the physical Q normalization.
"""
import sp_environment
import sys
sys.path.insert(0, str(sp_environment.WORK/'cluster-next-v2.nosync/work/so8-protected'))
import hw_scan as base
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache


@lru_cache(maxsize=16384)
def word_action(word, row, column):
    assert row != column and 0 <= row < 3 and 0 <= column < 3
    substitutions = {column: (1, row), 3+row: (-1, 3+column)}
    paths, traces = word
    result = defaultdict(int)
    for kind, sequences in enumerate((paths, traces)):
        for position, sequence in enumerate(sequences):
            for index, letter in enumerate(sequence):
                if letter not in substitutions:
                    continue
                coefficient, replacement = substitutions[letter]
                changed = sequence[:index]+(replacement,)+sequence[index+1:]
                sequences1 = sequences[:position]+(changed,)+sequences[position+1:]
                sign, target = base.canon(sequences1, traces) if kind == 0 else base.canon(paths, sequences1)
                if sign:
                    result[target] += sign*coefficient
    return {w:c for w,c in result.items() if c}


def action(expr, row, column, prime=None):
    result = defaultdict(int)
    for word, coefficient in expr.items():
        for target, factor in word_action(word, row, column).items():
            result[target] += coefficient*factor
    if prime is None:
        return {w:c for w,c in result.items() if c}
    return {w:int(c)%prime for w,c in result.items() if c%prime}


def root_factor(expr, row, column, weight, prime=None, shift=None):
    assert row < column
    shift = column-row if shift is None else shift
    result = dict(expr)
    raised = expr
    coefficient = 1
    for power in range(1, 1000):
        raised = action(raised, row, column, prime)
        if not raised:
            return result
        denominator = power*(weight+shift+power)
        assert denominator and (prime is None or denominator%prime)
        coefficient *= -Fraction(1, denominator) if prime is None else -pow(denominator, -1, prime)
        if prime is not None:
            coefficient %= prime
        lowered = raised
        for _ in range(power):
            lowered = action(lowered, column, row, prime)
        for word, value in lowered.items():
            result[word] = result.get(word, 0)+coefficient*value
        if prime is not None:
            result = {w:int(c)%prime for w,c in result.items() if c%prime}
        else:
            result = {w:c for w,c in result.items() if c}
    raise AssertionError('Nonterminating root string')


def project(expr, sector, prime=None, wrong_middle_shift=False):
    assert sector[1] >= sector[2] >= sector[3]
    for row, column in ((1,2), (0,2), (0,1)):
        expr = root_factor(expr, row, column, sector[row+1]-sector[column+1], prime,
                           1 if wrong_middle_shift and (row,column)==(0,2) else None)
    return expr


def insert(expr, pivots, prime):
    value = {w:int(c)%prime for w,c in expr.items() if c%prime}
    while value:
        word = min(value)
        coefficient = value[word]
        if word not in pivots:
            inverse = pow(coefficient, -1, prime)
            pivots[word] = {w:c*inverse%prime for w,c in value.items()}
            return True
        for w,c in pivots[word].items():
            value[w] = (value.get(w,0)-coefficient*c)%prime
            if not value[w]:
                del value[w]
    return False


def basis(sector, parity, prime):
    pivots = {}
    result = []
    for word in base.sector_basis(sector, parity):
        value = project({word:1}, sector, prime)
        if insert(value, pivots, prime):
            assert not action(value,0,1,prime) and not action(value,1,2,prime)
            result.append(value)
    return result
