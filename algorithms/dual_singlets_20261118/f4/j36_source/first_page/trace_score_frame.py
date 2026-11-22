"""Contract the corrected intermediate functional with two defining traces of the glueball."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
from itertools import combinations
import json
from math import lcm
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from flag_boundaries import HERE


def contractions(prepared):
    raw = (prepared/"queries.in").read_text().splitlines()
    skip = 1+int(raw[0])
    keys = [tuple(map(int,line.split())) for line in raw[skip+1:]]
    both = json.loads((prepared/"weights.json").read_bytes())
    coordinates = json.loads((prepared/"coordinate-rows.json").read_bytes())
    exact = json.loads((HERE.parent/"four-source-pairings.json").read_bytes())
    assert int(raw[skip]) == len(keys) == len(both) == 2985
    homogeneous = [Fraction(0)]*len(keys)
    for coefficient,row in zip(exact["etaK_weights"],coordinates,strict=True):
        for j,v in row:
            homogeneous[j] += Fraction(coefficient)*Fraction(v)
    weights = {}
    for k,(a,b,p,q),h in zip(keys,both,homogeneous,strict=True):
        for prime,value in ((1000003,p),(1000033,q)):
            assert h.numerator*pow(h.denominator,-1,prime)%prime == value
        weight = -4*Fraction(a)+2*Fraction(b)+h
        if weight:
            weights[k] = weight
    den = lcm(*(v.denominator for v in weights.values()))
    old = {k:int(v*den) for k,v in weights.items()}
    colour = json.loads(gzip.decompress((HERE/"inputs.json.gz").read_bytes()))
    tr = {(i,j):int(6*Fraction(v)) for i,j,v in colour["inverse"]}
    assert all(tr[i,j]==6*Fraction(v) for i,j,v in colour["inverse"])
    def insert_glue(row,n):
        out = Counter()
        for k,x in row.items():
            for i,u in enumerate(k[:n]):
                for j,v in enumerate(k[n:2*n]):
                    z = tr.get((u,v),0)
                    if z:
                        other = k[:i]+k[i+1:n+j]+k[n+j+1:]
                        out[other] += (-1)**(i+n+j-1)*x*z
        return {k:v for k,v in out.items() if v}
    single = insert_glue(old,5)
    twice = insert_glue(single,4)
    direct = Counter()
    for k,x in old.items():
        for i,j in combinations(range(5),2):
            for a,b in combinations(range(5,10),2):
                # The exterior square contributes minus twice the 2-by-2 minor.
                det = tr.get((k[i],k[a]),0)*tr.get((k[j],k[b]),0)-tr.get((k[i],k[b]),0)*tr.get((k[j],k[a]),0)
                if det:
                    other = tuple(t for r,t in enumerate(k) if r not in (i,j,a,b))
                    direct[other] -= 2*(-1)**(i+j+a+b-6)*det*x
    assert {k:v for k,v in direct.items() if v}==twice
    return dict(raw_coordinates=len(old),after_one_trace=len(single),after_two_traces=len(twice),
                denominator=den,entries=[[list(k),str(v)] for k,v in sorted(twice.items())])


if __name__ == "__main__":
    import sys
    sys.path.insert(0,str(HERE.parent))
    from bank_endpoints import prepare
    from verify import ARCHIVE
    from hashlib import sha256
    parser = argparse.ArgumentParser(description=__doc__)
    how = parser.add_mutually_exclusive_group(required=True)
    how.add_argument("--prepared",type=Path)
    how.add_argument("--archive",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    with TemporaryDirectory() as tmp:
        folder = args.prepared or Path(tmp)
        if args.archive:
            assert sha256(args.archive.read_bytes()).hexdigest()==ARCHIVE
            with ZipFile(args.archive) as z:
                prepare(z,folder)
        result = contractions(folder)
    data = json.dumps(result,separators=(",",":")).encode()
    args.out.write_bytes(gzip.compress(data,mtime=0))
    print({k:v for k,v in result.items() if k!="entries"})
