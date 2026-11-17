"""Partition an unchanged evaluation frame by points across native threads.

The threads use disjoint frame slices and thread-local native error storage.
No points, exterior components, coefficients or source recipes are changed.
"""
import sp_environment
import copy
import ctypes
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import native_su_derivation_v1 as derivation
import native_su_values_v1 as values
import native_su_profile_derivation_v1 as profiles

HERE = Path(__file__).resolve().parent
ENABLED = set()


def enable(prime):
    if prime in ENABLED:
        return
    reference = derivation.library(prime)
    library = ctypes.CDLL(str(HERE / ('su_derivation_d3_p%d_threadsafe_v1.dylib' % prime)))
    library.bmn_derivation_error.restype = ctypes.c_char_p
    library.bmn_derivation_evaluate.argtypes = reference.bmn_derivation_evaluate.argtypes
    library.bmn_derivation_evaluate.restype = ctypes.c_int
    library.trace_frame_error.restype = ctypes.c_char_p
    library.trace_frame_evaluate.argtypes = [ctypes.c_int] * 5 + [
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_int64), ctypes.c_int, ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_ubyte), ctypes.POINTER(ctypes.c_int64)]
    library.trace_frame_evaluate.restype = ctypes.c_int
    derivation.LIBRARIES[prime] = library
    ENABLED.add(prime)


def slices(frame, workers):
    key = (workers, frame.P)
    cache = getattr(frame, '_native_su_parallel_slices', {})
    if key not in cache:
        assert len(frame.degrees) == 1 and workers >= 1
        bounds = np.linspace(0, frame.P, min(workers, frame.P) + 1, dtype=int)
        result = []
        for start, stop in zip(bounds[:-1], bounds[1:]):
            item = copy.copy(frame)
            for name in list(item.__dict__):
                if name.startswith('_native_su_'):
                    delattr(item, name)
            item.P = int(stop - start)
            item.ncols = item.P * item.cpp
            item.L = {letter: (parity, data[start:stop]) for letter, (parity, data) in frame.L.items()}
            result.append(item)
        cache[key] = result
        frame._native_su_parallel_slices = cache
    return cache[key]


def apply(expressions, frame, evaluator, workers=8, *args):
    enable(frame.p)
    parts = slices(frame, workers)
    if len(parts) == 1:
        return evaluator(expressions, parts[0], *args)
    with ThreadPoolExecutor(max_workers=len(parts)) as pool:
        futures = [pool.submit(evaluator, expressions, part, *args) for part in parts]
        return np.concatenate([future.result() for future in futures], axis=1)


def evaluate(expressions, frame, workers=8):
    return apply(expressions, frame, values.evaluate, workers)


def evaluate_derivation(expressions, frame, workers=8):
    return apply(expressions, frame, derivation.evaluate, workers)


def evaluate_profile(expressions, frame, source_profile, target_profile, workers=8):
    return apply(expressions, frame, profiles.evaluate, workers, source_profile, target_profile)
