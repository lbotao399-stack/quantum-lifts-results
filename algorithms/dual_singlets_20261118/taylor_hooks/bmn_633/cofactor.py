"""Select a rational SO(7) cycle in the kernel of both Q0 and Cartan restriction."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from flint import nmod_mat
import numpy as np


PRIME = 2147483647
HERE = Path(__file__).resolve().parent


def read_matrix(path):
    with path.open("rb") as stream:
        rows, cols = map(int, np.fromfile(stream, dtype=np.uint32, count=2))
        matrix = np.fromfile(stream,dtype=np.uint32)
    assert matrix.size==rows*cols and int(matrix.max())<PRIME
    return matrix.reshape(rows,cols)


def pair(row, matrix):
    """Multiply over the prime field with bounded signed 64-bit dot products."""
    base = 1<<15
    hi, lo = divmod(np.asarray(row,dtype=np.int64),base)
    high, low = divmod(np.asarray(matrix,dtype=np.int64),base)
    assert len(row)<=4000 and len(row)*((PRIME//base)+1)**2 < 2**62
    return ((lo@low)%PRIME + ((hi@low+lo@high)%PRIME)*base
            + ((hi@high)%PRIME)*base**2)%PRIME


def column(values):
    return nmod_mat([[int(x)] for x in values],PRIME)


def entries(vector):
    return np.asarray([int(vector[i,0]) for i in range(vector.nrows())],dtype=np.int64)


def save_selected(path, rows, columns, folder):
    probes = [json.loads(line) for line in (folder/"probes.jsonl").read_bytes().splitlines()]
    generators = [json.loads(line) for line in (folder/"generators.jsonl").read_bytes().splitlines()]
    assert len(generators)==1837 and len(probes)>1837
    data = dict(prime=PRIME, ordered_rows=[probes[i] for i in rows],
                ordered_columns=[generators[j] for j in columns])
    raw = (json.dumps(data,separators=(",",":"))+"\n").encode()
    path.write_bytes(gzip.compress(raw,mtime=0))
    return hashlib.sha256(raw).hexdigest()


def run(folder, archive):
    maps = json.loads((HERE/"residual_so.json").read_bytes())
    assert maps["seed"]==317083 and maps["prime"]==PRIME
    b = read_matrix(folder/"incoming.bin")
    e = read_matrix(folder/"cochains.bin")
    a = read_matrix(folder/"closure_Cartan.bin")
    assert b.shape==(1837,1863) and e.shape==(1837,1837) and a.shape[1]==1837
    full = json.loads((HERE/"resolved.json").read_bytes())
    assert (full["SO7_incoming"],full["full_Cartan_character"]["highest_weight_dimension"],
            full["classical"]["SO7"]["exact_ranks"]["10"])==(1151,20,665)
    br, bc = (maps["incoming"][k] for k in ("rows","columns"))
    ar, ac = (maps["outgoing_and_Cartan"][k] for k in ("rows","columns"))
    assert len(br)==len(bc)==1151 and len(ar)==len(ac)==685
    bs = nmod_mat(b[np.ix_(br,bc)].tolist(),PRIME)
    small = nmod_mat(a[np.ix_(ar,ac)].tolist(),PRIME)
    assert int(bs.det())==maps["incoming"]["det"] and int(small.det())==maps["outgoing_and_Cartan"]["det"]

    # Project a deterministic functional off every incoming column.
    trial = np.asarray([pow(65537,i,PRIME) for i in range(len(e))],dtype=np.int64)
    dual = trial.copy()
    dual[br] = (dual[br]-entries(bs.transpose().solve(column(pair(trial,b[:,bc])))))%PRIME
    assert not pair(dual,b).any()
    on_generators = pair(dual,e)
    multipliers = entries(small.transpose().solve(column(on_generators[ac])))
    separations = (on_generators-pair(multipliers,a[ar]))%PRIME
    free, = np.nonzero(separations)
    assert len(free) and not separations[ac].any()
    j = int(free[0])
    cycle = np.zeros(e.shape[1],dtype=np.int64)
    cycle[j] = 1
    cycle[ac] = -entries(small.solve(column(a[ar,j])))%PRIME
    assert not pair(cycle,a.T).any()
    evaluation = pair(cycle,e.T)
    pairing = int(pair(dual,evaluation.reshape(-1,1))[0])
    assert pairing==int(separations[j])!=0
    digest = save_selected(archive,ar,ac+[j],folder)
    return dict(gauge="SO(7)",weight=[6,3,3],length=10,prime=PRIME,seed=317083,
                exact_incoming_rank=len(br),exact_joint_rank=len(ar),joint_minor_det=int(small.det()),
                outgoing_rows=sum(i<772 for i in ar),Cartan_rows=sum(i>=772 for i in ar),
                free_generator=j,free_generators_with_nonzero_pairing=len(free),
                nonzero_reduced_coefficients=int(np.count_nonzero(cycle)),integral_cofactor_pairing=pairing*int(small.det())%PRIME,
                free_coordinate_one_pairing=pairing,modular_dual=list(map(int,dual)),
                modular_cycle=[[int(i),int(cycle[i])] for i in np.flatnonzero(cycle)],
                exact_recipe_uncompressed_sha256=digest)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frames",type=Path,required=True)
    p.add_argument("--recipe",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__: p.error("enable assertions for verification")
    result=run(a.frames,a.recipe)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result["free_generator"],result["integral_cofactor_pairing"],flush=True)
