"""Three consecutive gauge-singlet cochain dimensions at every allowed source."""
from functools import cache
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from b2_schur_complex import multiplicity
from f4_singlet_characters import decomposition

@cache
def trace(parts):
    slots=[[],[]]
    for k in reversed(parts):
        min(slots,key=sum).append(k)
    left,right=[decomposition('F4',tuple(a)) for a in slots]
    return sum(v*right.get(w,0) for w,v in left.items())

def rows(p):
    for r in range(1,p-2):
        shape=(p-r,p-r)+(3,)*r
        degree=sum(shape)
        dims=[multiplicity(shape,n,trace) for n in (degree-2,degree-1,degree)]
        yield dict(power=p,order=r,shape=shape,length=degree-1,cochain_dimensions=dims)

if __name__=='__main__':
    print(json.dumps([z for p in (6,7) for z in rows(p)], indent=2))
