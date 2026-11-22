"""Verify all physical source grades, integer product minors, and the threshold return."""
import argparse
from collections import Counter
from functools import cache
from hashlib import sha256
import gzip
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import descendant
from scan import run


HERE=Path(__file__).resolve().parent
PRIMES=(1000003,1000033)
COUNTS={8:(0,0,0,0,3,0,0),9:(2,1,1,58,0,6,0)}


def classical(repo,ranks,rebuild,workers):
    native=repo/"algorithms/dual_singlets_20261118/taylor_hooks/native"
    for file,digest in json.loads((native/"sha256.json").read_bytes()).items():
        assert sha256((native/file).read_bytes()).hexdigest()==digest
    with TemporaryDirectory() as temp:
        scratch=Path(temp)
        if rebuild:
            engine=scratch/"engine"
            subprocess.run(["c++","-O3","-std=c++17","-pthread",native/"engine.cpp","-o",engine],check=True)
        for first,expected in COUNTS.items():
            for prime in PRIMES:
                if rebuild:
                    import contextlib
                    with contextlib.redirect_stdout(sys.stderr):
                        run(engine,scratch,first,range(7),prime,workers,900)
                for t in range(7):
                    key=f"Y{first}-t{t}-p{prime}"
                    row=ranks[key]
                    if rebuild:
                        fresh=json.loads((scratch/(key+".json")).read_bytes())
                        assert row==fresh["blocks"][0]
                    assert row["Y"]==first-t
                    assert row["F"]==first%2
                    assert row["dim"]-row["rank_in"]-row["rank_out"]==row["h"]==expected[t]
                    assert row==ranks[f"Y{first}-t{t}-p{PRIMES[0]}"]


def products(repo):
    sys.path.insert(0,str(repo/"algorithms/dual_singlets_20261118"))
    from flint import fmpz_mat
    from sp3_omega_product_replay import trace_factor,ordered,odd_jet
    from sp3_omega_products import raising,weight
    import mixed_incoming as m
    @cache
    def diagonal(spec):
        factors=trace_factor(spec,unitary=True)
        varied=Counter()
        diag=Counter()
        assert sum(spec[:5])==2
        for (word,),c in factors.items():
            for end,n in m.fw.Q_multi((word,)).items():
                varied[end]+=c*n
            jets=tuple((p,q,sum(1<<i for i in f)) for p,q,f in word)
            sign,labels=ordered(jets,odd_jet)
            if sign: diag[labels]+=2*sign*c
        assert not any(varied.values())
        assert all(n.denominator==1 for n in diag.values())
        return {w:int(n) for w,n in diag.items() if n}

    @cache
    def multiply(specs):
        if not specs: return {():1}
        out=Counter()
        for u,a in multiply(specs[:-1]).items():
            for v,b in diagonal(specs[-1]).items():
                sign,w=ordered(u+v,odd_jet)
                if sign: out[w]+=sign*a*b
        return {k:v for k,v in out.items() if v}

    records=[]
    for first,t in ((8,4),(9,1),(9,3),(9,5)):
        data=json.loads(gzip.decompress((HERE/f"products-{t}.json.gz").read_bytes()))
        wanted=(t,t,6-t,6-t,6-t)
        assert tuple(data["weight"])==wanted and data["length"]==first-t
        columns=[tuple(tuple(x) for x in slot) for (slot,) in data["columns"]]
        assert len(set(columns))==len(columns)==COUNTS[first][t]
        matrix=[]
        for recipe in data["recipes"]:
            raised=Counter(); row=Counter()
            for raw,c in recipe:
                specs=tuple(tuple(g) for g in raw)
                assert len(specs)*2==first-t and tuple(map(sum,zip(*(weight(g) for g in specs))))==wanted
                for key,n in raising(specs).items(): raised[key]+=c*n
                for w,n in multiply(specs).items(): row[w]+=c*n
            assert not any(raised.values())
            matrix.append([row[w] for w in columns])
        assert len(matrix)==len(columns)
        det=int(fmpz_mat(matrix).det())
        assert det and str(det)==data["minor"]
        records.append(dict(first_length=first,t=t,products=len(matrix),determinant=str(det)))
    return records


def verify(repo,rebuild,workers):
    ranks=json.loads((HERE/"ranks.json").read_bytes())
    classical(repo,ranks,rebuild,workers)
    product_results=products(repo)
    result=descendant.calculate()
    assert result["total_return"]==24883200
    result["single_jet_bracket_tests"]=descendant.verify_brackets()
    assert result==json.loads((HERE/"descendant.json").read_bytes())
    return dict(classical_source_dimensions={k:list(v) for k,v in COUNTS.items()},
                first_image_rank_in_KKLP_source=1,product_minors=product_results,
                returned_threshold=result,closed_source_dimensions_after_first=[3,67])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--rebuild",action="store_true")
    p.add_argument("--workers",type=int,default=2)
    a=p.parse_args()
    if not __debug__ or a.workers<1: p.error("enable assertions and use positive workers")
    print(json.dumps(verify(a.repo.resolve(),a.rebuild,a.workers),indent=2))
