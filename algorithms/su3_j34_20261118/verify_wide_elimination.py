"""Direct unblocked FLINT comparisons, degenerate panels, primal and dual witnesses."""
import json
import sys
import numpy as np
import wide_profile10 as w
sys.path.insert(0,str(w.f.v.maps.RUNTIME))

def run():
    rng=np.random.default_rng(917132601);P=w.P;records=[]
    for n,m,r,panel in [(9,12,0,3),(19,22,7,3),(9,7,7,2),(23,32,11,4),(18,19,17,31)]:
        A=rng.integers(P,size=(n,r),dtype=np.int64)@rng.integers(P,size=(r,m),dtype=np.int64)%P
        if n-r>=4:A[:3]=0;A[-1]=A[-2]
        assert int(w.f.matrix(A).rank())==r
        base=rng.integers(P,size=n)@A%P
        for i in range(2):
            b=base.copy();b[0]=(b[0]+i)%P; rec,cert=w.certify(w.fs.packed(A),w.fs.packed(b),panel)
            rank=int(w.f.matrix(np.vstack((A,b))).rank()); assert (rec['rank'],rec['augmented_rank'])==(r,rank)
            I,J=rec['pivot_rows'],rec['pivot_columns'];d=int(w.f.matrix(A[np.ix_(I,J)]).det()) if r else 1;assert d==rec['minor_determinant']
            if r<rank:
                cols=J+[rec['extra_column']];dd=int(w.f.matrix(np.vstack((A[I],b))[:,cols]).det());assert dd==rec['augmented_minor_determinant'] and int(b@cert%P)==rec['pairing'] and not(A@cert%P).any()
            else:np.testing.assert_array_equal(cert@A%P,b)
            records.append([n,m,panel,r,rank,rec['certificate']])
    return records
if __name__=='__main__':print(json.dumps(run()))
