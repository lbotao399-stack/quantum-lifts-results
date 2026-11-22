"""Check the rational saturation and which symmetric-group types the two slices detect."""
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import prod
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[2]))
from b2_schur_complex import centralizer, character, partitions


def multiplicity(shape, ordinary, signed):
    total=Fraction()
    for choice in product(*(partitions(k) for k in ordinary+signed)):
        sign=(-1)**sum(sum(c)-len(c) for c in choice[len(ordinary):])
        cycles=tuple(sorted(k for c in choice for k in c))
        total+=Fraction(sign*character(shape,cycles),prod(centralizer(c) for c in choice))
    assert total.denominator==1 and total>=0
    return int(total)


def verify():
    receipt=json.loads((HERE/"polarization_a8.json").read_bytes())
    assert receipt["accounting"].split("|")[1]=="COMPLETED" and receipt["accounting"].endswith("|0:0")
    for name,h in receipt["outputs"].items():
        assert sha256((HERE/name).read_bytes()).hexdigest()==h
    answers=[]
    for degrees in ((1,1,1),(2,1,1),(2,2,1)):
        n=sum(degrees)
        rows=[]
        for shape in sorted((p[::-1] for p in partitions(n)),reverse=True):
            a=multiplicity(shape,degrees,())
            b=multiplicity(shape,(1,),(n-1,))
            assert a+b>0
            rows.append([shape,a,b])
        for path in ("mixed_"+"".join(map(str,degrees)),f"dotted_{n-1}1"):
            d=json.loads((HERE/(path+".json")).read_bytes())
            for t in d["blocks"]:
                assert t["classical_rank_lower"]+t["known_cycle_rank_lower"]==t["columns"]
                assert t["unresolved_upper"]==0
            for key,value in d["total"].items():
                assert value==sum(t["orbit"]*t[key] for t in d["blocks"])
        answers.append(dict(degree=n,ordinary=degrees,one_scalar_dotted=n-1,types=rows))
    return answers


if __name__=="__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(verify(),indent=2))
