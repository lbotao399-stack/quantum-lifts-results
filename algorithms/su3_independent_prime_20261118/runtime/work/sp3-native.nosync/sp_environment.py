"""Resident dependencies and the existing SO/Sp trace-reversal quotient."""
import pathlib
import sys
from collections import defaultdict

HERE = pathlib.Path(__file__).resolve().parent
WORK = HERE.parent
for path in (WORK/'resident-deps.nosync',
             WORK/'cluster-next-v2.nosync/work/so8-bmn/vendor',
             WORK/'cluster-next-v2.nosync/work/upstream'):
    sys.path.insert(0, str(path))

from fortuity.full import words as fw
from fortuity.full.letters import is_ferm


def canonical(expr):
    # Same graded reversal as quantum_j18.so_canon. For symplectic matrices
    # transpose conjugation by the supplied form cancels around each trace.
    result = defaultdict(int)
    for multi, coefficient in expr.items():
        traces = []
        for word in multi:
            sign, value = fw.canon_word(word)
            reverse_sign, reverse = fw.canon_word(word[::-1])
            odd_count = sum(is_ferm(letter) for letter in word)
            reverse_sign *= (-1)**(len(word)+odd_count*(odd_count-1)//2)
            if not sign or (value == reverse and sign != reverse_sign):
                coefficient = 0
                break
            if reverse < value:
                value, sign = reverse, reverse_sign
            coefficient *= sign
            traces.append(value)
        if coefficient:
            sign, value = fw.canon_multi(traces)
            if sign:
                result[value] += coefficient*sign
    return {word: coefficient for word, coefficient in result.items() if coefficient}
