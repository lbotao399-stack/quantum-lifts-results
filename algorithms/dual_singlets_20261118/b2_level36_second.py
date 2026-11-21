"""Compute the second-order Taylor target of the fourteen-letter USp(4) BMN sources."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from b2_schur_complex import multiplicity
from taylor_hook_cohomology import saved, weights
from trace_descent import dimensions


PRIME = 2147483647
SHAPE = (5,5,1,1,1,1,1,1)
SECTORS = weights(5,6,True)
PREVIOUS = Path(__file__).resolve().parent/"taylor_hooks/cartan_products/symplectic_classical"


def dim(weight,length):
    return sum(dimensions((*weight[3:],*weight[:3]),length,4,"Sp").values())


def validate(data,length,i):
    w = SECTORS[i]
    sizes = [dim(w,length),dim(w,length+1)]
    assert (data["mode"],data["gauge"],data["sector"],data["length"],data["prime"]) == (
            "complete_classical_image","Sp2",list(w),length,PRIME)
    assert data["seed"] == 419903+1009*i+7919*length and data["dimensions"] == sizes
    minor = data["minor"]
    rank = minor["rank"]
    assert 0 <= rank <= min(sizes) and 0 < minor["det"] < PRIME
    for label,size in zip(("columns","rows"),sizes,strict=True):
        entries=minor[label]
        assert len(entries)==len(set(entries))==rank and all(0<=j<size for j in entries)
    return rank


def generate(engine,folder,tasks,threads):
    folder.mkdir(exist_ok=True,parents=True)
    env = dict(os.environ,PAR_THREADS=str(threads))
    for length,i in tasks:
        target = folder/f"image_{length}_{i}.json"
        if target.is_file():
            record = json.loads(target.read_bytes())
        else:
            seed = 419903+1009*i+7919*length
            command = [str(engine.resolve()),"Sp2",*map(str,SECTORS[i]),str(length),str(PRIME),str(seed)]
            record=json.loads(subprocess.check_output(command,env=env))
        rank=validate(record,length,i)
        saved(target,record)
        print(length,i,record["dimensions"],rank,flush=True)


def resolve(folder):
    images = {}
    for i,w in enumerate(SECTORS):
        file = next(PREVIOUS.glob(f"*/Sp2_5_6_{i}.json"))
        prior = json.loads(file.read_bytes())
        assert prior["sector"]==list(w) and prior["prime"]==PRIME and prior["blocks"][0]["Y"]==12
        images[12,i]=prior["blocks"][0]["rank_out"]
        for length in (13,14,15):
            data=json.loads((folder/f"image_{length}_{i}.json").read_bytes())
            images[length,i]=validate(data,length,i)
        images[16,i]=0
    blocks=[]
    for y in (13,14,15,16):
        cochains=multiplicity(SHAPE,y)
        assert sum((-1)**i*dim(w,y) for i,w in enumerate(SECTORS))==cochains
        ranks=[sum((-1)**i*images[k,i] for i in range(len(SECTORS))) for k in (y-1,y)]
        quotient=cochains-sum(ranks)
        assert min(*ranks,quotient)>=0
        blocks.append(dict(target_length=y,source_length=y+1,cochains=cochains,
                           incoming=ranks[0],outgoing=ranks[1],rational_cohomology_upper=quotient))
    return dict(group="USp(4)",source_weight=[10,4,4],order=2,partition=SHAPE,prime=PRIME,blocks=blocks)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder",type=Path,required=True)
    p.add_argument("--engine",type=Path)
    p.add_argument("--map",type=int,nargs=2,action="append",metavar=("LENGTH","INDEX"))
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    if not __debug__ or not (a.map or a.out) or bool(a.map)!=bool(a.engine):
        p.error("supply --map with --engine or --out with all certified maps")
    if a.map:
        assert all(y in (13,14,15) and i in range(len(SECTORS)) for y,i in a.map)
        generate(a.engine,a.folder,a.map,a.threads)
    if a.out:
        saved(a.out,resolve(a.folder))
