"""Inspect the trace-count filtration of an already certified defect."""
from inverse_support import *
from flint import nmod_mat
from fractions import Fraction
import math

start=time.monotonic()
folder=BASE_PACKAGE/'symplectic/sp5' if BASE_PACKAGE else HERE.parent/'fortuity_frontier_delivery_v2_20261123/symplectic/sp5'
data=json.loads(gzip.decompress((folder/'sp5_new_fortuitous_operator.json.gz').read_bytes()))
polynomial={pf.m.decode(w):Fraction(c) for c,w in data['terms']}
denominator=math.lcm(*(c.denominator for c in polynomial.values()))
defect=pf.m.q({w:int(c*denominator) for w,c in polynomial.items()})
sources=[]
for i in range(188):
    _,image=pickle.loads((pc.HERE/'sp5_q13_identity_rows'/f'row_{i}.pkl').read_bytes())
    sources.append(image)
records=[]
for cutoff in range(1,8):
    echelon=pf.Echelon();selected=[]
    for i,image in enumerate(sources):
        if echelon.add({w:c for w,c in image.items() if len(w)<=cutoff}):selected.append(i)
    rank=len(selected)
    independent=echelon.add({w:c for w,c in defect.items() if len(w)<=cutoff})
    record=dict(max_trace_factors=cutoff,boundary_rank=rank,augmented_rank=rank+int(independent),
                detects_class=independent,prime=1000003,seconds=time.monotonic()-start)
    records.append(record);print(json.dumps(record),flush=True)
(HERE/'trace_filtration_sp5.json').write_text(json.dumps(dict(records=records),indent=2)+'\n')
