"""Compare filtered elimination with the direct rank of every truncated matrix."""
from itertools import combinations_with_replacement
import numpy as np
from flint import nmod_mat

from filtered_rank import filtered_ranks


def verify():
    rng=np.random.default_rng(187311)
    matrices=intervals=0
    for prime in (3,101,1000003):
        for grades in range(1,7):
            for _ in range(125):
                rows=list(map(int,rng.integers(0,7,size=grades)))
                cols=list(map(int,rng.integers(0,7,size=grades)))
                blocks={}
                for a,b in combinations_with_replacement(range(grades),2):
                    if rng.integers(4):
                        rank=int(rng.integers(min(rows[a],cols[b])+1))
                        left=rng.integers(0,prime,(rows[a],rank))
                        right=rng.integers(0,prime,(rank,cols[b]))
                        blocks[a,b]=(left@right)%prime
                measured=filtered_ranks(blocks,rows,cols,prime)
                direct={}
                for lo in reversed(range(grades)):
                    for hi in range(lo,grades):
                        r=np.cumsum([0]+rows[lo:hi+1])
                        c=np.cumsum([0]+cols[lo:hi+1])
                        matrix=np.zeros((r[-1],c[-1]),dtype=np.int64)
                        for (a,b),block in blocks.items():
                            if lo<=a<=b<=hi:
                                matrix[r[a-lo]:r[a-lo+1],c[b-lo]:c[b-lo+1]]=block
                        direct[lo,hi]=nmod_mat(matrix.tolist(),prime).rank() if matrix.size else 0
                for a,b,value,count in measured:
                    assert value==direct[a,b]
                    assert count==value-direct.get((a+1,b),0)-direct.get((a,b-1),0)+direct.get((a+1,b-1),0)
                matrices+=1
                intervals+=len(measured)
    return matrices,intervals


if __name__=="__main__":
    assert __debug__
    print(*verify())
