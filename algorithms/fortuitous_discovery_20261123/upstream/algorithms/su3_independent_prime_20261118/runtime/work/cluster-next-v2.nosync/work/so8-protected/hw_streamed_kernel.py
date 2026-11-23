"""Exact formal kernels from streamed sparse coefficient constraints.

Every output coefficient is imposed. No sampled constraints or numerical
rank cutoff is used. Large incidence components avoid a dense matrix with
one row for every output word.
"""
from hw_blocked import *
import hw_blocked as blocked
from array import array

DENSE_ORIGINAL = blocked.DENSE_KERNEL


def as_array(M):
    return np.array([[int(M[i,j]) for j in range(M.ncols())]
                     for i in range(M.nrows())], dtype=np.int64).reshape(M.nrows(), M.ncols())


def bounded_kernel(images, p, chunk_rows=256, dense_entries=20000000, verbose=False):
    n = len(images)
    words = sorted({w for e in images for w,c in e.items() if int(c)%p})
    if n * len(words) <= dense_entries:
        return DENSE_ORIGINAL(images,p)
    start = time.monotonic()
    index = {w:i for i,w in enumerate(words)}
    rr, cc, vv = array('i'), array('i'), array('q')
    for j,e in enumerate(images):
        for w,c in e.items():
            if int(c)%p:
                rr.append(index[w]); cc.append(j); vv.append(int(c)%p)
    A = csr_matrix((np.frombuffer(vv,dtype=np.int64),
                    (np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),
                   shape=(len(words),n))
    assert n*(p-1)**2 < 2**63
    del index,rr,cc,vv
    H = np.eye(n,dtype=np.int64)
    changes = 0
    for offset in range(0,len(words),chunk_rows):
        values = np.asarray(A[offset:offset+chunk_rows] @ H.T) % p
        if not np.any(values):
            continue
        R, rank = nmod_mat(values.tolist(),p).rref()
        pivots = [next(j for j in range(R.ncols()) if R[i,j]) for i in range(rank)]
        pivset = set(pivots)
        free = [j for j in range(H.shape[0]) if j not in pivset]
        if not free:
            H = np.zeros((0,n),dtype=np.int64)
            break
        L = nmod_mat([[int(R[i,j]) for i in range(rank)] for j in free],p)
        correction = as_array(L * nmod_mat(H[pivots].tolist(),p))
        H = (H[free]-correction)%p
        changes += 1
        if verbose and (changes<=3 or changes%10==0):
            print(json.dumps(dict(stage='streamed_kernel', source_rows=n,
                                  output_words=len(words), imposed_through=offset+len(values),
                                  remaining_kernel=len(H), seconds=time.monotonic()-start)),flush=True)
    # Recheck all original coefficients against the final kernel.
    for offset in range(0,len(words),chunk_rows):
        assert not np.any(np.asarray(A[offset:offset+chunk_rows] @ H.T)%p)
    assert rank_mod(H,p)==len(H)
    if verbose:
        print(json.dumps(dict(stage='streamed_kernel_complete',source_rows=n,
                              output_words=len(words),kernel=len(H),rank=n-len(H),
                              seconds=time.monotonic()-start)),flush=True)
    return H,n-len(H)


def install():
    def kernel(images,p):
        return bounded_kernel(images,p,verbose=True)
    blocked.DENSE_KERNEL = kernel
