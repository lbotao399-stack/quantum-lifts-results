"""Expand both supercharges on every word of the raised SU(3) first equation."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
import gzip
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

try:
    from q1d_exact import Q1d_multi
    from fortuity.full.words import Q_multi
except ImportError:
    import product  # configure the archived trace runtime in the repository
    from q1d_exact import Q1d_multi
    from fortuity.full.words import Q_multi

LETTERS = [None]+[(a,b,t) for n in range(1,7) for a in range(n+1)
                 for b in range(n+1-a) for t in combinations(range(3),n-a-b)]
INDEX = {l:i for i,l in enumerate(LETTERS) if i}

@lru_cache(16384)
def pack_trace(t):
    return bytes(INDEX[l] for l in t)

def pack(w):
    return b"\0".join(pack_trace(t) for t in w)

def unpack(b):
    return tuple(tuple(LETTERS[i] for i in t) for t in b.split(b"\0"))

def export(path):
    from raised_first import transport
    x,y,r=transport()
    assert r["source_numerator_sha256"]=="8fcd4b0ebe148e20370cf452d24155bd464ed3ad4ff6712276c71addd119f875"
    assert r["correction_numerator_sha256"]=="62563243cb7a78faccc85c02a27b4efba8c71417fbe0d3c99bb8f221db274059"
    packet=[[[pack(w).hex(),v] for w,v in sorted(row.items())] for row in (x,y)]
    path.write_bytes(gzip.compress(json.dumps(packet,separators=(",",":")).encode(),mtime=0))

def apply(task):
    order,entries=task
    result=Counter()
    for h,v in entries:
        w=unpack(bytes.fromhex(h))
        image=Q1d_multi(w,3) if order else Q_multi(w)
        for t,c in image.items():
            result[pack(t)]+=v*c
    return {w:v for w,v in result.items() if v}

def verify(path,workers):
    raw=gzip.decompress(path.read_bytes())
    x,y=json.loads(raw)
    assert len(x)==34492 and len(y)==70141
    totals=[]
    with ProcessPoolExecutor(workers) as pool:
        for order,entries in ((1,x),(0,y)):
            row=Counter()
            tasks=((order,entries[i:i+128]) for i in range(0,len(entries),128))
            for partial in pool.map(apply,tasks):
                row.update(partial)
            totals.append({w:v for w,v in row.items() if v})
    rhs,lhs=totals
    assert rhs and lhs=={w:-v for w,v in rhs.items()}
    h=sha256()
    for w,v in sorted(rhs.items()):
        h.update(w.hex().encode()+b" "+str(Fraction(v)).encode()+b"\n")
    return dict(input_uncompressed_sha256=sha256(raw).hexdigest(),source_terms=len(x),
                first_correction_terms=len(y),nonzero_terms_on_each_side=len(rhs),
                residual_terms=0,Q1_source_sha256=h.hexdigest())

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("packet",type=Path)
    p.add_argument("--export",action="store_true")
    p.add_argument("--workers",type=int,default=1)
    args=p.parse_args()
    if not __debug__ or args.workers<1:
        p.error("enable assertions and use a positive worker count")
    if args.export:
        export(args.packet)
    else:
        print(json.dumps(verify(args.packet,args.workers),indent=2))
