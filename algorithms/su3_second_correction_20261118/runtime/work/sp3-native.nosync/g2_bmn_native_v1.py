"""G2 BMN frames and exact modular trace/Q0 evaluation in seven dimensions."""
from collections import defaultdict
from pathlib import Path
import ctypes
import numpy as np
from scipy.sparse import csr_matrix

import sp_census
from fortuity.hwcoh.evalg import Frame as OriginalFrame
from fortuity.bmn.words import QLET, contents_of_sector
from g2_full_classical_screen_v1 import g2_basis
from lie_root_dimension_v1 import multiplicity_parts

HERE = Path(__file__).resolve().parent
LIBRARIES = {}


class Frame(OriginalFrame):
    def __init__(self, degree, points, seed, prime, k=None):
        super().__init__(7, points, degree if k is None else k, [degree], seed=seed, p=prime, group='SU')
        self.group = 'G2'
        basis, _ = g2_basis()
        rng = np.random.default_rng(seed)
        for letter in sorted(self.L):
            odd = int(3 <= letter <= 5)
            coefficients = rng.integers(0, prime, (points, self.k if odd else 1, 14), dtype=np.int64)
            self.L[letter] = odd, np.einsum('...t,tij->...ij', coefficients, basis) % prime

    def values(self, words):
        assert all(not paths for paths, traces in words)
        return self.eval_words([traces for paths, traces in words], dtype=np.int64)


def configure(prime):
    sp_census.configure(3, prime)
    def expected(sector, parity):
        assert parity == 0
        answer = defaultdict(int)
        for content in contents_of_sector(sector):
            answer[sum(content[3:6])] += multiplicity_parts(
                tuple(sorted(content[:3] + (content[6],), reverse=True)),
                tuple(sorted(content[3:6], reverse=True)), 'G2')
        return dict(answer)
    sp_census.base.expected = expected
    sp_census.base.Frame = Frame
    sp_census.base.exact_hw.cache_clear()


def library(prime):
    if prime not in LIBRARIES:
        native = ctypes.CDLL(str(HERE / ('bmn_derivation_rank_v1_d7_p%d.dylib' % prime)))
        ints = ctypes.POINTER(ctypes.c_int)
        offsets = ctypes.POINTER(ctypes.c_uint64)
        values = ctypes.POINTER(ctypes.c_int64)
        encoded = ctypes.POINTER(ctypes.c_ubyte)
        native.trace_frame_evaluate.argtypes = [ctypes.c_int] * 5 + [ints, offsets, values, ctypes.c_int, offsets, encoded, values]
        native.trace_frame_evaluate.restype = ctypes.c_int
        native.trace_frame_error.restype = ctypes.c_char_p
        native.bmn_derivation_evaluate.argtypes = [ctypes.c_int] * 5 + [ints, offsets, values, ctypes.c_int, ints, ctypes.c_int, offsets, encoded, values]
        native.bmn_derivation_evaluate.restype = ctypes.c_int
        native.bmn_derivation_error.restype = ctypes.c_char_p
        LIBRARIES[prime] = native
    return LIBRARIES[prime]


def word_values(words, frame, derivative=False):
    assert frame.group == 'G2' and len(frame.degrees) == 1
    if not words:
        return np.zeros((0, frame.ncols), dtype=np.int64)
    degree = frame.degrees[0]
    cached = getattr(frame, '_g2_native_configuration', None)
    if cached is None:
        ids = {letter: index + 1 for index, letter in enumerate(sorted(frame.L))}
        parities = np.zeros(len(ids) + 1, dtype=np.int32)
        offsets = np.zeros(len(ids) + 1, dtype=np.uint64)
        matrices, offset = [], 0
        for letter, index in ids.items():
            odd, values = frame.L[letter]
            assert values.shape[-2:] == (7, 7)
            parities[index], offsets[index] = odd, offset
            flat = np.ascontiguousarray(values, dtype=np.int64).reshape(-1)
            matrices.append(flat)
            offset += len(flat)
        rules = np.asarray([(ids[source], coefficient, ids[replacement[0]], ids[replacement[1]])
                            for source, replacements in sorted(QLET.items())
                            for coefficient, replacement in replacements], dtype=np.int32)
        cached = ids, parities, offsets, np.concatenate(matrices), rules
        frame._g2_native_configuration = cached
    ids, parities, offsets, matrices, rules = cached
    encoded, word_offsets = bytearray(), [0]
    for word in words:
        count = sum(3 <= letter <= 5 for trace in word for letter in trace)
        assert abs(count - degree) == 1 if derivative else count == degree
        for trace in word:
            encoded.extend(ids[letter] for letter in trace)
            encoded.append(0)
        word_offsets.append(len(encoded))
    encoded = np.frombuffer(bytes(encoded), dtype=np.uint8)
    word_offsets = np.asarray(word_offsets, dtype=np.uint64)
    output = np.empty((len(words), frame.ncols), dtype=np.int64)
    native = library(frame.p)
    args = [frame.k, frame.P, degree, frame.p, len(ids),
            parities.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            matrices.ctypes.data_as(ctypes.POINTER(ctypes.c_int64))]
    if derivative:
        args += [len(rules), rules.ctypes.data_as(ctypes.POINTER(ctypes.c_int))]
    args += [len(words), word_offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
             encoded.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
             output.ctypes.data_as(ctypes.POINTER(ctypes.c_int64))]
    status = (native.bmn_derivation_evaluate if derivative else native.trace_frame_evaluate)(*args)
    assert status == 0, (native.bmn_derivation_error if derivative else native.trace_frame_error)().decode()
    return output


def evaluate(expressions, frame, derivative=False):
    degree = frame.degrees[0]
    def keep(word):
        count = sum(3 <= letter <= 5 for trace in word[1] for letter in trace)
        return abs(count - degree) == 1 if derivative else count == degree
    words = sorted({word for expression in expressions for word in expression if keep(word)})
    assert all(not word[0] for word in words)
    if not words:
        return np.zeros((len(expressions), frame.ncols), dtype=np.int64)
    indices = {word: index for index, word in enumerate(words)}
    rows, columns, coefficients = [], [], []
    for row, expression in enumerate(expressions):
        for word, coefficient in expression.items():
            if word in indices:
                rows.append(row)
                columns.append(indices[word])
                coefficients.append(int(coefficient) % frame.p)
    matrix = csr_matrix((np.asarray(coefficients, dtype=np.int64), (rows, columns)),
                        shape=(len(expressions), len(words)))
    assert int(max(np.diff(matrix.indptr), default=0)) * (frame.p - 1) ** 2 < 2 ** 63
    return np.asarray(matrix @ word_values([word[1] for word in words], frame, derivative)) % frame.p
