"""Pinned BMN extremal projection through the full-letter native projector.

The BMN middle dual letter is minus theta_0 theta_2. Conjugating by
that existing change of basis reproduces hw_extremal exactly. The native
library owns mutable work buffers: use separate processes, never threads.
"""
from pathlib import Path
import ctypes
import numpy as np

HERE = Path(__file__).resolve().parent
LIBRARY = None


def library():
    global LIBRARY
    if LIBRARY is None:
        value = ctypes.CDLL(str(HERE/'bmn_projection_v1.dylib'))
        value.full_projection.argtypes = [ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.c_int,
            ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_ubyte),
            ctypes.POINTER(ctypes.c_int64)]
        value.full_projection.restype = ctypes.c_int
        value.full_projection_error.restype = ctypes.c_char_p
        value.full_projection_size.restype = ctypes.c_uint64
        value.full_projection_copy.argtypes = [ctypes.POINTER(ctypes.c_ubyte),
            ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int64)]
        LIBRARY = value
    return LIBRARY


def project(expr, sector, prime=1000003):
    assert prime in (1000003, 1000033)
    assert len(sector) == 4 and sector[1] >= sector[2] >= sector[3]
    if not expr:
        return {}
    data = np.asarray([(0, 0, mask) for mask in (1, 2, 4, 6, 5, 3, 7)], dtype=np.int32)
    raw, offsets, coefficients = bytearray(), [0], []
    for (paths, traces), coefficient in expr.items():
        assert not paths
        middle_duals = 0
        for trace in traces:
            assert all(0 <= letter < 7 for letter in trace)
            raw.extend(letter+1 for letter in trace)
            raw.append(0)
            middle_duals += trace.count(4)
        offsets.append(len(raw))
        coefficients.append(int(coefficient)*(-1 if middle_duals % 2 else 1) % prime)
    raw = np.frombuffer(bytes(raw), dtype=np.uint8)
    offsets = np.asarray(offsets, dtype=np.uint64)
    coefficients = np.asarray(coefficients, dtype=np.int64)
    full_sector = np.asarray((0, 0, 2*sector[1], 2*sector[2], 2*sector[3]), dtype=np.int32)
    native = library()
    count = native.full_projection(prime, 7,
        data.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        full_sector.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), len(coefficients),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        raw.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    assert count >= 0, native.full_projection_error().decode()
    raw = np.empty(native.full_projection_size(), dtype=np.uint8)
    offsets = np.empty(count+1, dtype=np.uint64)
    coefficients = np.empty(count, dtype=np.int64)
    native.full_projection_copy(raw.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    raw = raw.tobytes()
    result = {}
    for index, coefficient in enumerate(coefficients):
        key = raw[int(offsets[index]):int(offsets[index+1])]
        traces = tuple(tuple(letter-1 for letter in trace) for trace in key[:-1].split(b'\0'))
        middle_duals = sum(trace.count(4) for trace in traces)
        result[((), traces)] = int(coefficient)*(-1 if middle_duals % 2 else 1) % prime
    return result
