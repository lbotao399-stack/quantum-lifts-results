"""Evaluate the pinned BMN QLET rules by the graded matrix product rule."""
import sp_environment
from fortuity.bmn.words import QLET
from pathlib import Path
import ctypes
import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
LIBRARIES = {}


def library(prime):
    if prime not in LIBRARIES:
        value = ctypes.CDLL(str(HERE/('bmn_derivation_rank_v1_d6_p%d.dylib' % prime)))
        value.bmn_derivation_error.restype = ctypes.c_char_p
        value.bmn_derivation_evaluate.argtypes = [ctypes.c_int]*5 + [
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_uint64),
            ctypes.POINTER(ctypes.c_int64), ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.c_int, ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_ubyte),
            ctypes.POINTER(ctypes.c_int64)]
        value.bmn_derivation_evaluate.restype = ctypes.c_int
        LIBRARIES[prime] = value
    return LIBRARIES[prime]


def word_images(words, frame):
    assert frame.group == 'Sp' and frame.N == 3 and len(frame.degrees) == 1
    prime, degree = frame.p, frame.degrees[0]
    assert prime in (1000003, 1000033)
    if not words:
        return np.zeros((0, frame.ncols), dtype=np.int64)
    cached = getattr(frame, '_native_bmn_derivation_configuration', None)
    if cached is None:
        alphabet = sorted(frame.L)
        ids = {letter: index+1 for index, letter in enumerate(alphabet)}
        parities = np.zeros(len(alphabet)+1, dtype=np.int32)
        offsets = np.zeros(len(alphabet)+1, dtype=np.uint64)
        matrices, offset = [], 0
        for letter, index in ids.items():
            parity, values = frame.L[letter]
            assert values.shape[-2:] == (6, 6) and values.shape[0] == frame.P
            parities[index], offsets[index] = parity, offset
            values = np.ascontiguousarray(values, dtype=np.int64).reshape(-1)
            matrices.append(values)
            offset += len(values)
        rules = []
        for source, replacements in sorted(QLET.items()):
            for coefficient, replacement in replacements:
                assert len(replacement) == 2
                rules.append((ids[source], coefficient, ids[replacement[0]], ids[replacement[1]]))
        rules = np.asarray(rules, dtype=np.int32).reshape(-1, 4)
        cached = ids, parities, offsets, np.concatenate(matrices), rules
        frame._native_bmn_derivation_configuration = cached
    ids, parities, offsets, matrices, rules = cached
    encoded, word_offsets = bytearray(), [0]
    for word in words:
        assert abs(sum(3 <= letter <= 5 for trace in word for letter in trace)-degree) == 1
        for trace in word:
            encoded.extend(ids[letter] for letter in trace)
            encoded.append(0)
        word_offsets.append(len(encoded))
    encoded = np.frombuffer(bytes(encoded), dtype=np.uint8)
    word_offsets = np.asarray(word_offsets, dtype=np.uint64)
    output = np.empty((len(words), frame.ncols), dtype=np.int64)
    native = library(prime)
    status = native.bmn_derivation_evaluate(
        frame.k, frame.P, degree, prime, len(ids),
        parities.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        matrices.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)), len(rules),
        rules.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(words),
        word_offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        encoded.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        output.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    assert status == 0, native.bmn_derivation_error().decode()
    return output


def evaluate(expressions, frame):
    degree = frame.degrees[0]
    words = sorted({traces for expr in expressions for paths, traces in expr
                    if not paths and abs(sum(3 <= letter <= 5 for trace in traces for letter in trace)-degree) == 1})
    assert all(not paths for expr in expressions for paths, traces in expr)
    if not words:
        return np.zeros((len(expressions), frame.ncols), dtype=np.int64)
    indices = {word: index for index, word in enumerate(words)}
    rows, columns, coefficients = [], [], []
    for index, expr in enumerate(expressions):
        for (_, word), coefficient in expr.items():
            if word in indices:
                rows.append(index)
                columns.append(indices[word])
                coefficients.append(int(coefficient) % frame.p)
    assert max(map(len, expressions), default=0)*(frame.p-1)**2 < 2**63
    matrix = csr_matrix((np.asarray(coefficients, dtype=np.int64), (rows, columns)),
                        shape=(len(expressions), len(words)))
    return np.asarray(matrix @ word_images(words, frame)) % frame.p
