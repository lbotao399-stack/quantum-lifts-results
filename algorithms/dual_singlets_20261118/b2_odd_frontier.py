"""Bound the first odd-length BMN generator using the published SO(5) ranks."""
import argparse
import hashlib
import json
from pathlib import Path

import mixed_incoming as m
from trace_descent import dimensions


def audit(path, out):
    raw=path.read_bytes()
    public=json.loads(raw)
    assert public["jmax"]==28 and not public["checks_failed"]
    rows={(tuple(d["lam"]),d["twoJ"]):d for d in public["columns"]}
    assert len(rows)==len(public["columns"])==public["sectors"]
    witness=[]
    for total in range(1,15):
        for x in range(1,total+1):
            for y in range(min(x,total-x)+1):
                z=total-x-y
                if not 0<=z<=y:
                    continue
                weight=(x,y,z)
                for length in range(x,total+1):
                    if length%2==0:
                        continue
                    spin=total-length
                    dim=sum(dimensions((0,0,*weight),length,4,"Sp").values())
                    sector=rows.get((weight,spin))
                    if not dim:
                        assert sector is None or sector["m"]==0
                        continue
                    assert sector["check"] and sector["j"]==2*total and sector["m"]==dim
                    preceding=sector["r_up"]
                    next_record=rows.get((weight,spin-1))
                    if next_record is None:
                        assert not sum(dimensions((0,0,*weight),length+1,4,"Sp").values())
                    else:
                        assert next_record["check"]
                    following=0 if next_record is None else next_record["r_up"]
                    assert dim-preceding-following==sector["h"]==0
                    witness.append([list(weight),length,dim,preceding,following])
    assert witness
    data=dict(public_table_sha256=hashlib.sha256(raw).hexdigest(),
              public_rows=len(rows),nonempty_odd_blocks=len(witness),
              rational_odd_cohomology_through_level_28=0,nonempty_blocks=witness)
    m.write(out,data)
    print("Odd-length blocks with adjacent ranks exhausting the rational character",len(witness))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    audit(args.table,args.out)
