"""Integer exterior-algebra witness at bidegree (2,2), in the scalar trace normalization."""
from itertools import permutations
import json
from scalar_upper_bound import compute as scalar


def wedge(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():
            if i&j:continue
            sign=(-1)**sum((j % (1<<k)).bit_count() for k in range(4) if i&(1<<k))
            m=i|j;out[m]=out.get(m,0)+sign*x*y
    return {i:x for i,x in out.items() if x}


def compute():
    out={}
    for g,control in scalar().items():
        gram=control['trace_pairings'];bilinear={}
        for i in range(2):
            for j in range(2):
                if gram[i][j]:bilinear[(1<<i)|(4<<j)]=gram[i][j]
        square=wedge(bilinear,bilinear);assert square=={15:-8}
        det=control['determinant'];coeff=det*square[15]
        assert coeff=={'SO(7)':64,'USp(6)':-64}[g]
        out[g]=dict(scalar_trace_determinant=det,gaugino_trace_by_exterior_mask=bilinear,
                    squared_gaugino_trace_by_exterior_mask=square,witness_top_coefficient=coeff)
    return out

if __name__=='__main__':print(json.dumps(compute(),indent=2))
