"""Exact native evaluation of the existing Sp(2) matrix/exterior frames."""
import sp_environment
from pathlib import Path
import ctypes
import numpy as np
from scipy.sparse import csr_matrix

HERE=Path(__file__).resolve().parent
_library=ctypes.CDLL(str(HERE/'trace_frame_native.dylib'))
_library.trace_frame_error.restype=ctypes.c_char_p
_library.trace_frame_evaluate.argtypes=[ctypes.c_int]*5+[
    ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_uint64),ctypes.POINTER(ctypes.c_int64),ctypes.c_int,
    ctypes.POINTER(ctypes.c_uint64),ctypes.POINTER(ctypes.c_ubyte),ctypes.POINTER(ctypes.c_int64)]
_library.trace_frame_evaluate.restype=ctypes.c_int


def evaluate_words(words,frame):
    assert len(frame.degrees)==1 and frame.p==1000003
    if not words:return np.zeros((0,frame.ncols),dtype=np.int64)
    cached=getattr(frame,'_native_trace_configuration',None)
    if cached is None:
        alphabet=sorted(frame.L);ids={L:i+1 for i,L in enumerate(alphabet)}
        parities=np.zeros(len(alphabet)+1,dtype=np.int32);offsets=np.zeros(len(alphabet)+1,dtype=np.uint64)
        matrices=[];offset=0
        for i,L in enumerate(alphabet,1):
            degree,values=frame.L[L];assert values.shape[-2:]==(4,4) and values.shape[0]==frame.P
            parities[i]=degree;offsets[i]=offset;values=np.ascontiguousarray(values,dtype=np.int64).reshape(-1)
            matrices.append(values);offset+=len(values)
        cached=ids,parities,offsets,np.concatenate(matrices);frame._native_trace_configuration=cached
    ids,parities,offsets,matrices=cached;raw=bytearray();word_offsets=[0]
    for word in words:
        for trace in word:raw.extend(ids[L] for L in trace);raw.append(0)
        word_offsets.append(len(raw))
    raw=np.frombuffer(bytes(raw),dtype=np.uint8);word_offsets=np.asarray(word_offsets,dtype=np.uint64)
    values=np.empty((len(words),frame.ncols),dtype=np.int64)
    status=_library.trace_frame_evaluate(frame.k,frame.P,frame.degrees[0],frame.p,len(ids),
        parities.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
        matrices.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),len(words),
        word_offsets.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),raw.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte)),
        values.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)))
    assert status==0,_library.trace_frame_error().decode()
    return values


def evaluate(expressions,frame):
    words=sorted({word for expr in expressions for word in expr})
    if not words:return np.zeros((len(expressions),frame.ncols),dtype=np.int64)
    indices={word:i for i,word in enumerate(words)};rows=[];columns=[];coefficients=[]
    for i,expr in enumerate(expressions):
        for word,c in expr.items():rows.append(i);columns.append(indices[word]);coefficients.append(int(c)%frame.p)
    assert max(map(len,expressions),default=0)*(frame.p-1)**2<2**63
    matrix=csr_matrix((np.array(coefficients,dtype=np.int64),(rows,columns)),shape=(len(expressions),len(words)))
    return np.asarray(matrix@evaluate_words(words,frame))%frame.p
