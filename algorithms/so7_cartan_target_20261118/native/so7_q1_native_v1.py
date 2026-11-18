"""Native SO7 Q1 with the existing rational derivative pair rule."""
import sp_environment
from pathlib import Path
from fractions import Fraction
import ctypes
import numpy as np
from q1d_sp_exact import pair_terms
from fortuity.full.letters import is_ferm

LIB = ctypes.CDLL(str(Path(__file__).with_suffix('.dylib')))
I = ctypes.POINTER(ctypes.c_int)
U64 = ctypes.POINTER(ctypes.c_uint64)
I64 = ctypes.POINTER(ctypes.c_int64)
U8 = ctypes.POINTER(ctypes.c_ubyte)
LIB.so7_q1_error.restype = ctypes.c_char_p
LIB.so7_q1_configure.argtypes = [ctypes.c_int, ctypes.c_int, I, ctypes.c_int, I]
LIB.so7_q1_configure.restype = ctypes.c_int
LIB.so7_q1_apply.argtypes = [ctypes.c_int, U64, U8, I64]
LIB.so7_q1_apply.restype = ctypes.c_int
LIB.so7_q1_size.restype = ctypes.c_uint64
LIB.so7_q1_copy.argtypes = [U8, U64, I64]
CONFIG = None

def configure(initial, prime=1000003):
    global CONFIG
    initial = tuple(sorted(set(initial)))
    if CONFIG is not None and CONFIG['initial'] == initial and CONFIG['prime'] == prime:
        return CONFIG
    alphabet = set(initial); rules = []
    for left in initial:
        for right in initial:
            for coefficient, first, second in pair_terms(left, right):
                coefficient = Fraction(coefficient)
                value = coefficient.numerator * pow(coefficient.denominator, -1, prime) % prime
                if value:
                    alphabet.update((first, second)); rules.append((left, right, value, first, second))
    alphabet = sorted(alphabet); assert len(alphabet) <= 255
    ids = {letter:i+1 for i,letter in enumerate(alphabet)}
    parities = np.asarray([int(is_ferm(letter)) for letter in alphabet], dtype=np.int32)
    raw = np.asarray([(ids[l], ids[r], c, ids[a], ids[b]) for l,r,c,a,b in rules], dtype=np.int32).reshape((-1,5))
    code = LIB.so7_q1_configure(prime, len(alphabet), parities.ctypes.data_as(I), len(raw), raw.ctypes.data_as(I))
    assert code == 0, LIB.so7_q1_error().decode()
    CONFIG = dict(initial=initial, prime=prime, alphabet=alphabet, ids=ids, rules=len(raw))
    return CONFIG

def apply(expression, prime=1000003, initial=None):
    if not expression:
        return {}
    if initial is None:
        initial = {letter for word in expression for trace in word for letter in trace}
    config = configure(initial, prime)
    raw = bytearray(); offsets=[0]; coefficients=[]
    for word, coefficient in expression.items():
        for trace in word:
            raw.extend(config['ids'][letter] for letter in trace);raw.append(0)
        offsets.append(len(raw))
        coefficient = Fraction(coefficient)
        coefficients.append(coefficient.numerator * pow(coefficient.denominator, -1, prime) % prime)
    raw = np.frombuffer(bytes(raw), dtype=np.uint8)
    offsets = np.asarray(offsets, dtype=np.uint64)
    coefficients = np.asarray(coefficients, dtype=np.int64)
    count = LIB.so7_q1_apply(len(coefficients), offsets.ctypes.data_as(U64), raw.ctypes.data_as(U8), coefficients.ctypes.data_as(I64))
    assert count >= 0, LIB.so7_q1_error().decode()
    output = np.empty(LIB.so7_q1_size(), dtype=np.uint8)
    offsets = np.empty(count+1, dtype=np.uint64); coefficients=np.empty(count, dtype=np.int64)
    LIB.so7_q1_copy(output.ctypes.data_as(U8), offsets.ctypes.data_as(U64), coefficients.ctypes.data_as(I64))
    data = output.tobytes(); alphabet=config['alphabet'];result={}
    for i,c in enumerate(coefficients):
        encoded=data[int(offsets[i]):int(offsets[i+1])]
        word=tuple(tuple(alphabet[j-1] for j in trace) for trace in encoded[:-1].split(b'\0'))
        result[word]=int(c)
    return result
