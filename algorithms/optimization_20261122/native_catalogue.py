"""Optional native front end for the unchanged F4 degree-fourteen catalogue."""
import ctypes
from functools import cache
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile

import numpy as np


@cache
def library():
    source = Path(__file__).with_name("dihedral_catalogue.cpp")
    digest = sha256(source.read_bytes()).hexdigest()[:16]
    directory = Path(tempfile.gettempdir()) / ("quantum-lifts-dihedral-" + digest)
    directory.mkdir(exist_ok=True)
    binary = directory / "dihedral.so"
    if not binary.exists():
        with tempfile.NamedTemporaryFile(dir=directory, suffix=".so", delete=False) as temporary:
            target = Path(temporary.name)
        try:
            subprocess.run(["c++", "-O3", "-std=c++17", "-shared", "-fPIC",
                            str(source), "-o", str(target)], check=True)
            target.replace(binary)
        finally:
            target.unlink(missing_ok=True)
    native = ctypes.CDLL(str(binary))
    native.dihedral_error.restype = ctypes.c_char_p
    native.dihedral_catalogue.argtypes = [ctypes.c_int,
                                         ctypes.POINTER(ctypes.c_int),
                                         ctypes.POINTER(ctypes.c_size_t)]
    native.dihedral_catalogue.restype = ctypes.POINTER(ctypes.c_uint64)
    native.dihedral_free.argtypes = [ctypes.c_void_p]
    native.dihedral_canonical14.argtypes = [ctypes.POINTER(ctypes.c_uint32),
                                          ctypes.POINTER(ctypes.c_uint32),
                                          ctypes.c_size_t]
    return native


def fixed_content_catalogue(counts):
    """Sorted base-four uint64 codes for all nonzero dihedral trace orbits."""
    counts = tuple(counts)
    if not 1 <= len(counts) <= 4 or any(type(c) is not int or c < 0 for c in counts):
        raise ValueError("one through four nonnegative integer counts required")
    if not 2 <= sum(counts) <= 31:
        raise ValueError("total trace length must be 2 through 31")
    content = (ctypes.c_int * len(counts))(*counts)
    size = ctypes.c_size_t()
    native = library()
    pointer = native.dihedral_catalogue(len(counts), content, ctypes.byref(size))
    if not pointer:
        raise RuntimeError(native.dihedral_error().decode())
    try:
        return np.ctypeslib.as_array(pointer, shape=(size.value,)).copy()
    finally:
        native.dihedral_free(pointer)


@cache
def catalogue():
    """Exact replacement for single_trace_flag.catalogue, including dtype."""
    return fixed_content_catalogue((4, 4, 4, 2)).astype(np.int32)


def canonical14(codes):
    """Exact replacement for single_trace_flag.canonical on base-four codes."""
    array = np.asarray(codes)
    if array.dtype.kind not in "iu" or np.any(array < 0) or np.any(array >= 1 << 28):
        raise ValueError("codes must be integers between zero and 4**14 - 1")
    flat = np.ascontiguousarray(array, dtype=np.uint32)
    output = np.empty_like(flat)
    library().dihedral_canonical14(flat.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                                  output.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
                                  flat.size)
    return output.astype(array.dtype).reshape(array.shape)
