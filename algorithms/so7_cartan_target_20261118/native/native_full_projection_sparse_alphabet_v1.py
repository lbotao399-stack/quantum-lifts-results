"""Native full projector using an alphabet closed under the actual letter actions."""
import sp_environment
from fortuity.full.letters import letters_up_to
from functools import lru_cache
from pathlib import Path
import ctypes
import numpy as np

HERE=Path(__file__).resolve().parent
_lib=ctypes.CDLL(str(HERE/'full_projection_native.dylib'))
_lib.full_projection.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.c_int,
    ctypes.POINTER(ctypes.c_uint64),ctypes.POINTER(ctypes.c_ubyte),ctypes.POINTER(ctypes.c_int64)]
_lib.full_projection.restype=ctypes.c_int
_lib.full_projection_error.restype=ctypes.c_char_p
_lib.full_projection_size.restype=ctypes.c_uint64
_lib.full_projection_copy.argtypes=[ctypes.POINTER(ctypes.c_ubyte),ctypes.POINTER(ctypes.c_uint64),ctypes.POINTER(ctypes.c_int64)]


@lru_cache(None)
def alphabet(total_derivatives):
    letters=sorted(L for L in letters_up_to(total_derivatives,total_derivatives) if L[0]+L[1]<=total_derivatives)
    assert len(letters)<=255
    data=np.asarray([(a,b,sum(1<<i for i in S)) for a,b,S in letters],dtype=np.int32)
    return letters,{L:i+1 for i,L in enumerate(letters)},data


def project(expr,sector,prime=1000003):
    assert prime==1000003
    if not expr:return {}
    derivatives=max(L[0]+L[1] for word in expr for trace in word for L in trace)
    letters,ids,data=alphabet(derivatives);raw=bytearray();offsets=[0];coefficients=[]
    for word,c in expr.items():
        for trace in word:raw.extend(ids[L] for L in trace);raw.append(0)
        offsets.append(len(raw));coefficients.append(int(c)%prime)
    raw=np.frombuffer(bytes(raw),dtype=np.uint8);offsets=np.asarray(offsets,dtype=np.uint64)
    coefficients=np.asarray(coefficients,dtype=np.int64);sector=np.asarray(sector,dtype=np.int32)
    count=_lib.full_projection(prime,len(letters),data.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        sector.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),len(coefficients),
        offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),raw.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    assert count>=0,_lib.full_projection_error().decode()
    size=_lib.full_projection_size();out=np.empty(size,dtype=np.uint8);offsets=np.empty(count+1,dtype=np.uint64);coefficients=np.empty(count,dtype=np.int64)
    _lib.full_projection_copy(out.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),coefficients.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    raw=out.tobytes();result={}
    for i,c in enumerate(coefficients):
        key=raw[int(offsets[i]):int(offsets[i+1])]
        word=tuple(tuple(letters[index-1] for index in trace) for trace in key[:-1].split(b'\0'))
        result[word]=int(c)
    return result
