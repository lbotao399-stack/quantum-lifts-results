"""Reconstruct the integer cofactor recipe at two primes from its saved trace words and probes."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from flint import nmod_mat
import numpy as np


def vector(items):
    return [len(items),*items]


def lists(items):
    return [len(items),*(i for row in items for i in vector(row))]


def encode(data):
    rows,cols=data["ordered_rows"],data["ordered_columns"]
    yield [len(rows),len(cols)]
    for g in cols:
        yield [g["skeleton"],*lists(g["cols"]),*lists(g["edges"]),*vector(g["mu"]),*lists(g["tr"])]
    for row in rows:
        yield [row["kind"],row["skeleton"],len(row["letters"])]
        for kind,flavour,plus,minus,matrices in row["letters"]:
            yield [kind,flavour,plus,minus,len(matrices)]
            yield from matrices


def replay(recipe, check, engine):
    body=gzip.decompress(recipe.read_bytes())
    report=json.loads(check.read_bytes())
    assert hashlib.sha256(body).hexdigest()==report["exact_recipe_uncompressed_sha256"]
    data=json.loads(body)
    supplied=("\n".join(" ".join(map(str,row)) for row in encode(data))+"\n").encode()
    out={}
    with tempfile.TemporaryDirectory() as td:
        dest=Path(td)/"rows.bin"
        for prime in (2147483647,1000033):
            subprocess.run([str(engine.resolve()),str(prime),str(dest)],input=supplied,check=True)
            values=np.fromfile(dest,dtype=np.uint32).reshape(685,686)
            assert values.max()<prime
            base=nmod_mat(values[:,:-1].tolist(),prime)
            determinant=int(base.det())
            assert determinant
            coeff=base.solve(nmod_mat(values[:,-1:].tolist(),prime))
            free_coordinate_one=[(-int(coeff[i,0]))%prime for i in range(685)]+[1]
            if prime==report["prime"]:
                minor=json.loads((check.parent/"residual_so.json").read_bytes())["outgoing_and_Cartan"]
                ordered=minor["columns"]+[report["free_generator"]]
                archived=dict(report["modular_cycle"])
                assert free_coordinate_one==[archived.get(i,0) for i in ordered]
                assert determinant==report["joint_minor_det"]
            out[str(prime)]=dict(determinant=determinant,nonzero_coefficients=sum(bool(v) for v in free_coordinate_one),
                                 modular_coefficients_sha256=hashlib.sha256(
                                     json.dumps(free_coordinate_one,separators=(",",":")).encode()).hexdigest())
    return dict(exact_recipe_uncompressed_sha256=report["exact_recipe_uncompressed_sha256"],
                ordered_integer_rows=dict(Counter(r["kind"] for r in data["ordered_rows"])),
                trace_generators=len(data["ordered_columns"]),primes=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("recipe","check","engine","out"):
        p.add_argument("--"+name,type=Path,required=True)
    a=p.parse_args()
    if not __debug__:p.error("enable assertions")
    a.out.write_text(json.dumps(replay(a.recipe,a.check,a.engine),indent=2)+"\n")
