"""Compare the old-kernel update with direct FLINT elimination, including zero growth."""
import json
import sys
from pathlib import Path

import numpy as np
import four_source_refit as r

sys.path.insert(0, str(r.v.maps.RUNTIME))


def controls():
    P = r.P
    packed = r.fs.packed
    rng = np.random.default_rng(32977513)
    pick = lambda shape: rng.integers(P, size=shape, dtype=np.int64)
    checks = []
    for old_rank, growth in ((1,0),(3,0),(2,1),(4,3),(5,5)):
        rows, width, extra = 13, 8, 7
        I = list(map(int, rng.permutation(rows)[:old_rank]))
        J = list(map(int, rng.permutation(width)[:old_rank]))
        U, V = pick((rows,old_rank)), pick((old_rank,width))
        U[I] = np.eye(old_rank,dtype=np.int64); V[:,J] = np.diag(np.arange(2,2+old_rank))
        A = U@V%P; c=pick(rows); y=c@A%P
        X = pick((rows,growth)); X[I]=0
        newrows=[i for i in reversed(range(rows)) if i not in I][:growth]
        X[newrows] = np.eye(growth,dtype=np.int64)
        Y = pick((growth,extra)); Y[:,:growth]=np.eye(growth,dtype=np.int64)
        H=(U@pick((old_rank,extra))+X@Y)%P
        good=(c@H + pick(growth)@Y)%P
        for consistent in (True,False):
            u=good.copy(); u[-1]=(int(u[-1])+int(not consistent))%P
            B=np.hstack((A,H)); b=np.r_[y,u]
            direct=[int(r.f.matrix(B).rank()),int(r.f.matrix(np.vstack((B,b))).rank())]
            d=int(r.f.matrix(A[np.ix_(I,J)]).det())
            info,kind,cert=r.extension(packed(A),packed(y),packed(H),packed(u),c,I,J,d)
            assert direct==[old_rank+growth,old_rank+growth+(not consistent)]==[info['rank'], info['augmented_rank']]
            rows2=I+info['extra_pivot_rows'];cols2=J+[width+j for j in info['extra_pivot_columns']]
            assert int(r.f.matrix(B[np.ix_(rows2,cols2)]).det())==info['full_minor_determinant']
            cert=cert.astype(np.int64)
            if consistent:
                assert kind=='coefficients'; np.testing.assert_array_equal(cert@B%P,b)
            else:
                assert kind=='separator' and int(b@cert%P)==info['pairing']!=0; assert not (B@cert%P).any()
            checks.append(dict(old_rank=old_rank,growth=growth,direct_ranks=direct,certificate=kind))
    return checks

if __name__=='__main__':
    print(json.dumps(controls(),indent=2))
