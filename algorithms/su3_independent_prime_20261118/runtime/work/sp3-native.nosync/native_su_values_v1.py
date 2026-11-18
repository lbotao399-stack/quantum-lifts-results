"""Cyclic SU(3) trace values using the already compiled native exterior engine."""
import sp_environment
from native_su_derivation_v1 import library
from fortuity.full.evalg import nferm
import ctypes
import numpy as np
from scipy.sparse import csr_matrix


def word_values(words, frame):
    assert frame.group == 'SU' and frame.N == 3 and len(frame.degrees) == 1
    degree = frame.degrees[0]
    if not words:
        return np.zeros((0, frame.ncols), dtype=np.int64)
    cached = getattr(frame, '_native_su_value_configuration', None)
    if cached is None:
        ids = {letter: index + 1 for index, letter in enumerate(sorted(frame.L))}
        parities = np.zeros(len(ids) + 1, dtype=np.int32)
        offsets = np.zeros(len(ids) + 1, dtype=np.uint64)
        matrices, offset = [], 0
        for letter, index in ids.items():
            parity, values = frame.L[letter]
            assert values.shape[-2:] == (3, 3) and values.shape[0] == frame.P
            parities[index], offsets[index] = parity, offset
            values = np.ascontiguousarray(values, dtype=np.int64).reshape(-1)
            matrices.append(values)
            offset += len(values)
        cached = ids, parities, offsets, np.concatenate(matrices)
        frame._native_su_value_configuration = cached
    ids, parities, offsets, matrices = cached
    encoded, word_offsets = bytearray(), [0]
    for word in words:
        assert nferm(word) == degree
        for trace in word:
            encoded.extend(ids[letter] for letter in trace)
            encoded.append(0)
        word_offsets.append(len(encoded))
    encoded = np.frombuffer(bytes(encoded), dtype=np.uint8)
    word_offsets = np.asarray(word_offsets, dtype=np.uint64)
    output = np.empty((len(words), frame.ncols), dtype=np.int64)
    native = library(frame.p)
    native.trace_frame_evaluate.argtypes = [ctypes.c_int] * 5 + [
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_int64), ctypes.c_int, ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_ubyte), ctypes.POINTER(ctypes.c_int64)]
    native.trace_frame_evaluate.restype = ctypes.c_int
    native.trace_frame_error.restype = ctypes.c_char_p
    status = native.trace_frame_evaluate(
        frame.k, frame.P, degree, frame.p, len(ids),
        parities.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        matrices.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)), len(words),
        word_offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        encoded.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        output.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    assert status == 0, native.trace_frame_error().decode()
    return output


def evaluate(expressions, frame):
    words = sorted({word for expr in expressions for word in expr if nferm(word) == frame.degrees[0]})
    if not words:
        return np.zeros((len(expressions), frame.ncols), dtype=np.int64)
    indices = {word: index for index, word in enumerate(words)}
    rows, columns, coefficients = [], [], []
    for index, expr in enumerate(expressions):
        for word, coefficient in expr.items():
            if word in indices:
                rows.append(index)
                columns.append(indices[word])
                coefficients.append(int(coefficient) % frame.p)
    assert max(map(len, expressions), default=0) * (frame.p - 1) ** 2 < 2 ** 63
    matrix = csr_matrix((np.asarray(coefficients, dtype=np.int64), (rows, columns)),
                        shape=(len(expressions), len(words)))
    return np.asarray(matrix @ word_values(words, frame)) % frame.p
