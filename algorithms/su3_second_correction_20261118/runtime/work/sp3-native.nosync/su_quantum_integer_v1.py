"""Exact rational Q1/Q2 in SU(N), with checked integer native accumulation."""
import sp_environment
from pathlib import Path
import sys
sys.path.insert(0, str(sp_environment.WORK / 'cluster-next-v2.nosync/work/so8-protected'))
from q1d_exact import pair_terms, Q1d_multi
from q2d_exact import triple_rule, jet_terms, Q2d_multi
from fortuity.full.letters import is_ferm
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import lcm
import ctypes
import struct

HERE = Path(__file__).resolve().parent
LIB = ctypes.CDLL(str(HERE / 'su_quantum_integer_v1.dylib'))
LIB.su_integer_error.restype = ctypes.c_char_p
LIB.su_integer_init.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
LIB.su_integer_rule.argtypes = [ctypes.c_int] * 4 + [ctypes.c_int64] + [ctypes.c_int] * 2
LIB.su_integer_apply.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_int64]
LIB.su_integer_finish.argtypes = [ctypes.POINTER(ctypes.c_uint64)]
LIB.su_integer_finish.restype = ctypes.POINTER(ctypes.c_ubyte)


def species_rules(key, order):
    result = defaultdict(Fraction)
    if order == 1:
        for coefficient, first, second in pair_terms(*key):
            result[first, second] += coefficient
    else:
        left, middle, right = key
        for coefficient, theta_first, theta_second in triple_rule(left[2], middle[2], right[2]):
            for weight, jet_first, jet_second in jet_terms(tuple(letter[:2] for letter in key)):
                result[jet_first + (theta_first,), jet_second + (theta_second,)] += Fraction(-coefficient, 6) * weight
    return {pair: value for pair, value in result.items() if value}


def act(expression, order, N=3):
    assert order in (1, 2)
    if not expression:
        return {}
    keys, alphabet = set(), set()
    for word in expression:
        letters = [letter for trace in word for letter in trace]
        alphabet.update(letters)
        if order == 1:
            keys.update((letters[i], letters[j]) for i, j in combinations(range(len(letters)), 2))
        else:
            for triple in combinations(range(len(letters)), 3):
                for middle in triple:
                    left, right = [i for i in triple if i != middle]
                    keys.add((letters[left], letters[middle], letters[right]))
    rules = {key: species_rules(key, order) for key in keys}
    for rule in rules.values():
        alphabet.update(letter for pair in rule for letter in pair)
    alphabet = sorted(alphabet)
    assert len(alphabet) <= 255
    ids = {letter: i + 1 for i, letter in enumerate(alphabet)}
    parities = (ctypes.c_int * (len(alphabet) + 1))(*([0] + [int(is_ferm(x)) for x in alphabet]))
    assert LIB.su_integer_init(N, len(alphabet), parities) == 0, LIB.su_integer_error()
    input_denominator = lcm(*(Fraction(c).denominator for c in expression.values()))
    rule_denominator = lcm(*(Fraction(c).denominator for rule in rules.values() for c in rule.values()))
    for key, rule in rules.items():
        positions = (ids[key[0]], 0, ids[key[1]]) if order == 1 else tuple(ids[x] for x in key)
        for (first, second), coefficient in rule.items():
            integer = coefficient * rule_denominator
            assert integer.denominator == 1 and abs(integer) < 2 ** 63
            assert LIB.su_integer_rule(order, *positions, int(integer), ids[first], ids[second]) == 0, LIB.su_integer_error()
    for word, coefficient in expression.items():
        integer = Fraction(coefficient) * input_denominator
        assert integer.denominator == 1 and abs(integer) < 2 ** 63
        encoded = bytes(i for trace in word for i in (*(ids[x] for x in trace), 0))
        assert LIB.su_integer_apply(order, encoded, len(encoded), int(integer)) == 0, LIB.su_integer_error()
    size = ctypes.c_uint64()
    pointer = LIB.su_integer_finish(ctypes.byref(size))
    data = ctypes.string_at(pointer, size.value)
    result, offset = {}, 0
    while offset < len(data):
        length, coefficient = struct.unpack_from('<Iq', data, offset)
        offset += 12
        encoded = data[offset:offset + length]
        offset += length
        assert encoded[-1] == 0 and coefficient
        word = tuple(tuple(alphabet[i - 1] for i in trace) for trace in encoded[:-1].split(b'\0'))
        assert word not in result
        result[word] = Fraction(coefficient, input_denominator * rule_denominator)
    assert offset == len(data)
    return result


def reference(expression, order, N=3):
    result = defaultdict(Fraction)
    operation = Q1d_multi if order == 1 else Q2d_multi
    for word, coefficient in expression.items():
        for target, value in operation(word, N, 'SU').items():
            result[target] += coefficient * value
    return {word: coefficient for word, coefficient in result.items() if coefficient}
