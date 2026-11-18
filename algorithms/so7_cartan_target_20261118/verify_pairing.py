"""Check the rank bounds, target covector and ordered determinant manifests."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from flint import nmod_mat
from check_incoming import check_target, cochain_dimension

HERE=Path(__file__).resolve().parent
P=1000003


def load(name, digest):
    raw=(HERE/name).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==digest
    return json.loads(gzip.decompress(raw))


def verify():
    target=check_target()
    with np.load(HERE/"target_frame.npz",allow_pickle=False) as z:
        rank=nmod_mat(z["boundaries"].tolist(),P).rank()
    assert rank==target["boundary_rows"]==1443 and target["SF_pairing"]==1
    successor=load("successor_minor.json.gz","ce1a6ef5bba22c0c19317170d0810d4b0c45ee52484b0f54624512b1279591e0")
    assert len(successor["words"])==len(set(successor["columns"]))==2307
    assert successor["determinant"]==488963
    incoming=load("incoming_minor.json.gz","6fbdb485472d3e9c6385e15a8d8bfc2cf690642a00f15dfb5e96321f44a4ed20")
    assert incoming["prime"]==P and len(incoming["panels"])==212
    columns=[]
    det=1
    for panel in incoming["panels"]:
        assert len(panel["words"])==len(panel["columns"]) and 0<panel["determinant"]<P
        columns.extend(panel["columns"])
        det=det*panel["determinant"]%P
    assert len(columns)==len(set(columns))==12526 and sum(incoming["widths"]) in columns
    assert det==incoming["determinant"]==942245
    dim=cochain_dimension((3,3,2,2,2),11)
    assert dim==14832 and dim-len(successor["words"])==len(columns)-1
    return dict(prime=P,target_boundary_rank=rank,target_pairing=1,
        successor_minor=[len(successor["words"]),successor["determinant"]],
        augmented_minor=[len(columns),det],classical_rank_over_Q=dim-len(successor["words"]),
        d1_rank_to_SF_over_Q=1)


if __name__=="__main__":
    print(json.dumps(verify(),indent=2))
