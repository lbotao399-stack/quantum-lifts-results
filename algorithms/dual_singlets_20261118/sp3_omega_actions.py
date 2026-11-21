"""Use the one-dimensional Sp(3) target to project the direct first and second supercharges."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import mixed_incoming as m
from dual_charge_bands import coordinates
from so8_level24_chain import solve_image
from trace_descent import action, basis, dimensions, evaluate, highest_weights, settings


def load_page(path):
    z = json.loads(gzip.decompress(path.read_bytes()))
    assert (z["weight"], z["length"], z["prime"], z["group"], z["colours"]) == (
        [3,3,5,0,0], 9, m.P, "Sp", 6)
    assert z["modular_cohomology"] == z["representative"]["pairing"] == 1
    return z


def centre_basis(page):
    rows = highest_weights(tuple(page["weight"]),page["length"])
    return [rows[int(degree)][i] for degree,frame in page["centre_frames"].items()
            for i in frame["minor_rows"]]


def signed_frame(weight, length, seed):
    rows = highest_weights(weight,length)
    dims = dimensions(weight,length,6,"Sp")
    alphabet = sorted({l for content in m.contents(weight,length) for l in content})
    _, frames = basis(rows,dims,"Sp",6,seed,alphabet)
    return frames


def rational_cycle(path, out):
    page = load_page(path)
    weight = tuple(page["weight"])
    chosen = centre_basis(page)
    residues = page["representative"]["cycle_basis_coefficients"]
    assert len(chosen) == len(residues)
    outgoing = page["outgoing"]["rows"]
    support = [i for i in outgoing if residues[i]]
    special, = set(i for i,c in enumerate(residues) if c).difference(outgoing)
    assert residues[special] == 1
    frames = signed_frame(weight,10,240191)
    images = [m.q(chosen[i]) for i in support+[special]]
    values = evaluate(images,frames,6,"Sp")
    print("rational solve",values.shape,flush=True)
    vector,info = solve_image(values[:-1],[-int(x) for x in values[-1]])
    coefficients = {special: Fraction(1), **{support[i]:c for i,c in vector.items()}}
    for i,c in coefficients.items():
        assert c.denominator % m.P and c.numerator*pow(c.denominator,-1,m.P)%m.P == residues[i]
    assert all(i in coefficients or c==0 for i,c in enumerate(residues))
    result = dict(page_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  cycle_basis_coefficients=[[i,str(c)] for i,c in sorted(coefficients.items())],
                  signed_frames=frames, exact_closed_coordinates=values.shape[1],
                  nonzero_residuals=0, modular_pairing=1, **info)
    m.write(out/"rational_cycle.json.gz",result)
    return result


def direct(path, order, out, workers):
    settings()
    page=load_page(path)
    grade=3-order
    weight=(grade,grade,8-grade,3-grade,3-grade)
    length=8+order
    src=highest_weights(weight,length)
    rows=[r for d in sorted(src) for r in src[d]]
    lib=m.library(out,6,m.P)
    target=(3,3,5,0,0)
    letters=sorted({l for content in m.contents(target,9) for l in content})
    assert set(page["centre_frames"])=={"4"}
    probes=page["centre_frames"]["4"]["samples"]
    functional=np.array(page["representative"]["target_functional"],dtype=np.int64)
    assert len(functional)*(m.P-1)**2 < 2**63
    result=np.zeros(len(rows),np.int64)
    batch_size=96
    checkpoint=out/f"direct{order}_checkpoint.json"
    start=0
    if checkpoint.is_file():
        cached=json.loads(checkpoint.read_bytes())
        assert cached["page_sha256"]==hashlib.sha256(path.read_bytes()).hexdigest()
        result[:len(cached["entries"])]=cached["entries"]
        start=len(cached["entries"])
    for lo in range(start,len(rows),batch_size):
        images=[action(r,order,6,"Sp") for r in rows[lo:lo+batch_size]]
        values=coordinates(images,4,letters,probes,6,"Sp",lib,workers)
        hi=lo+len(images)
        result[lo:hi]=values@functional%m.P
        if (lo//batch_size)%12==0 or hi==len(rows):
            m.write(checkpoint,dict(page_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                    entries=result[:hi].tolist()))
            print("Q",order,hi,len(rows),"nonzero",int(np.count_nonzero(result[:hi])),flush=True)
    record=dict(page_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),prime=m.P,
                order=order,source_weight=weight,source_length=length,formal_rows=len(rows),
                target_dimension=len(functional),entries=result.tolist(),nonzero=int(np.count_nonzero(result)))
    m.write(out/f"direct{order}.json.gz",record)
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--page",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--task",choices=("cycle","first","second"),required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    assert __debug__
    args.out.mkdir(exist_ok=True,parents=True)
    if args.task=="cycle":
        rational_cycle(args.page,args.out)
    else:
        direct(args.page,1 if args.task=="first" else 2,args.out,args.workers)
