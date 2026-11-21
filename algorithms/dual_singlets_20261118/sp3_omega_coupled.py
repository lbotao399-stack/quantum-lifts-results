"""Extend a dual representative of omega through the two incoming equations."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import mixed_incoming as m
from dual_charge_bands import coordinates
from sp3_omega_actions import load_page
from sp3_omega_top_dual import Echelon, all_rows
from trace_descent import action, dimensions, settings


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def samples(dims, seed):
    result = {}
    for degree, size in dims.items():
        points = math.ceil((size+12)/(degree+1))
        result[degree] = [(min(8,points-i), seed+317*degree+6151*i,
                           list(range(min(8,points-i)*(degree+1)))) for i in range(0,points,8)]
    return result


def frame_width(probes):
    return sum(len(cols) for blocks in probes.values() for _,_,cols in blocks)


def evaluate(rows, probes, alphabet, library, workers):
    return np.column_stack([coordinates(rows, degree, alphabet, lists, 6, "Sp", library, workers)
                            for degree,lists in probes.items()])


def profile(grade):
    weight = (grade,grade,8-grade,3-grade,3-grade)
    source_length = 11-grade
    alphabet = sorted({l for content in m.contents(weight,source_length+1) for l in content})
    return weight, source_length, alphabet


def direct(path, order, page_hash):
    data = read(path)
    weight, length, _ = profile(3-order)
    assert (data["order"], data["page_sha256"], data["prime"], data["source_weight"],
            data["source_length"]) == (order, page_hash, m.P, list(weight), length)
    assert len(data["entries"]) == data["formal_rows"] and sum(bool(v) for v in data["entries"]) == data["nonzero"]
    return np.asarray(data["entries"], dtype=np.int64)


def extend(grade, direct_entries, previous, library, workers, out):
    weight, length, alphabet = profile(grade)
    seeds = {2:932597, 1:402759}
    probes = samples(dimensions(weight,length+1,6,"Sp"), seeds[grade])
    width = frame_width(probes)
    assert width * (m.P-1)**2 < 2**63
    rows = all_rows(weight, length)
    assert len(rows) == len(direct_entries)
    # Keep the evaluation matrix in node memory. A whole-system equality can then
    # be verified after elimination, including rows where q0 vanishes identically.
    table = np.empty((len(rows),width),dtype=np.uint32)
    rhs = -direct_entries.copy() % m.P
    space = Echelon(width, rhs=True)
    if previous is not None:
        prev_probes, prev_vector = previous
        _, _, prev_alphabet = profile(grade+1)
    for lo in range(0,len(rows),192):
        batch = rows[lo:lo+192]
        hi = lo+len(batch)
        classical = evaluate([m.q(row) for row in batch], probes, alphabet, library, workers)
        table[lo:hi] = classical
        if previous is not None:
            first = evaluate([action(row,1,6,"Sp") for row in batch],
                             prev_probes, prev_alphabet, library, workers)
            rhs[lo:hi] = (rhs[lo:hi] - first@prev_vector) % m.P
        space.add(np.column_stack((classical,rhs[lo:hi])),list(range(lo,hi)))
        if (lo//192)%12==0 or hi==len(rows):
            msg=dict(grade=grade,processed=hi,source_rows=len(rows),partial_rank=space.rank,
                     target_coordinates=width,has_conflict=space.witness is not None)
            m.write(out/f"grade{grade}_progress.json",msg)
            print(msg,flush=True)
        if space.witness is not None:
            witness = np.zeros(hi,dtype=np.int64)
            for i,c in space.witness.items():
                witness[i]=c
            assert hi*(m.P-1)**2 < 2**63
            for j in range(0,width,256):
                assert not np.any(witness @ table[:hi,j:j+256].astype(np.int64) % m.P)
            assert int(witness@rhs[:hi] % m.P)==1
            return dict(grade=grade, solved=False, formal_rows=len(rows), tested_rows=hi,
                        target_coordinates=width, samples=probes,
                        conflict=[[i,c] for i,c in sorted(space.witness.items())]), None
    vector = space.solution()
    for lo in range(0,len(rows),256):
        end=min(lo+256,len(rows))
        assert np.array_equal(table[lo:end].astype(np.int64)@vector % m.P,rhs[lo:end])
    result=dict(grade=grade, solved=True, formal_rows=len(rows), verified_rows=len(rows),
                target_coordinates=width, rank=space.rank, nonzero_residuals=0,
                samples=probes, functional=vector.tolist(),
                rhs_nonzero=int(np.count_nonzero(rhs)),
                rhs_sha256=hashlib.sha256(rhs.astype("<u4").tobytes()).hexdigest())
    return result, (probes,vector)


def run(page_path, first_path, second_path, out, workers):
    settings()
    load_page(page_path)
    page_hash=hashlib.sha256(page_path.read_bytes()).hexdigest()
    paths=[first_path,second_path]
    out.mkdir(parents=True,exist_ok=True)
    library=m.library(out,6,m.P)
    previous=None
    for grade,path in zip((2,1),paths,strict=True):
        order=3-grade
        record,previous=extend(grade,direct(path,order,page_hash),previous,library,workers,out)
        record.update(prime=m.P,page_sha256=page_hash,direct_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        m.write(out/f"dual_grade{grade}.json.gz",record)
        if not record["solved"]:
            # A conflict for a fixed covector needs the full homogeneous dual freedom.
            return False
    return True


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for opt in ("page","first","second","out"):
        parser.add_argument("--"+opt,type=Path,required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    assert __debug__
    run(args.page,args.first,args.second,args.out,args.workers)
