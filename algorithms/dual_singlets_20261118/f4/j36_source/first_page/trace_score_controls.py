"""Compare native full-length forward action to the transposed, shortened trace calculation."""
import argparse
from fractions import Fraction
from functools import cache
import gzip
from itertools import combinations,product
import json
from math import factorial,prod
from pathlib import Path
import sys
from zipfile import ZipFile

from sympy import Matrix

from flag_boundaries import HERE
from trace_source_score import GLUE, PAIRS, native_data, section, templates

sys.path.insert(0,str(HERE.parent))
from first_action import FirstAction, inputs
from verify import ARCHIVE
from hashlib import sha256


def lifted(m):
    jets = [(0,0,1<<i) for i in range(3)]+[(1,0,0),(0,1,0)]+[(0,0,(1<<i)+(1<<j)) for i,j in PAIRS]
    out = []
    for term in product(*(section(f) if f!=GLUE else [(GLUE,1)] for f in m)):
        out.append(dict(profile="Psi",coefficient=str(prod(c for _,c in term)),
                        factors=[(1,[jets[x] for x in w]) for w,_ in term]))
    return out


def calculate(archive):
    assert sha256(archive.read_bytes()).hexdigest()==ARCHIVE
    with ZipFile(archive) as z:
        tables,_,raw = inputs(z)
    keys = sorted({tuple(key) for family in raw for key,_ in family})
    data = json.loads(gzip.decompress((HERE/"trace_source_cycles.json.gz").read_bytes()))
    tr,rep,omega = native_data()
    @cache
    def minor(p,m):
        return int(Matrix([[3*int(tr[a,b]) for b in m] for a in p]).det())
    @cache
    def numeric(w):
        a=rep[w[0]]
        for i in w[1:]:
            a=a@rep[i]
        return int(a.trace())
    def short(key,k,polys):
        n=5-k
        total = 0
        for i in combinations(range(5),k):
            for j in combinations(range(5,10),k):
                det = minor(tuple(key[t] for t in i),tuple(key[t] for t in j))
                if not det:
                    continue
                sign = (-1)**(sum(i+j)-k*(2*k-1)+k*(k-1)//2)
                plus = tuple(key[t] for t in range(5) if t not in i)
                minus = tuple(key[t] for t in range(5,10) if t not in j)
                for h,u in enumerate(plus):
                    for t,v in enumerate(minus):
                        emitted = (-1)**(h+n+t-1)
                        rest = plus[:h]+plus[h+1:]+minus[:t]+minus[t+1:]
                        for a,b,c in omega.get((u,v),()):
                            labels=(*key[10:],a,b,*rest)
                            for p,y in polys.items():
                                total += det*sign*emitted*c*y*prod(numeric(tuple(labels[s] for s in w)) for w in p)
        return Fraction(total*factorial(k),4*2**(13-2*k))

    checks = []
    key=(9,11,12,18,42,19,29,45,46,47,0,4,51)
    assert key in keys
    for k,i in ((2,147),(3,31),(4,0)):
        m=tuple(map(tuple,data["monomials"][i]))
        power,polys=templates(m)
        assert power==k
        value=short(key,k,polys)
        direct = FirstAction(tables,lifted(m)).coefficient(key)
        assert value and value==direct, (k,i,key,value,direct)
        checks.append(dict(traces_removed=k,source=i,raw_colour_query=key,
                           forward_coefficient=str(direct),transposed_coefficient=str(value)))
    return checks


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive",type=Path)
    args=parser.parse_args()
    if not __debug__:
        parser.error("enable assertions")
    print(json.dumps(calculate(args.archive),indent=2))
