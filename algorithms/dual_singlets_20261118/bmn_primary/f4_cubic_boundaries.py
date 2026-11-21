"""Insert the field-strength and scalar--fermion pair into the Albert cubic contractions."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mixed_incoming as m
from f4_bmn_double_jordan import inherited
from f4_bmn_primary import LETTERS, SIGNS
from f4_colour import frame
from f4_joined_boundaries import native
from f4_scalar_nine import act, alternate, on_traceless, pairing, tensors


def bracket(x, y, prime):
    return (x@y-y@x)%prime


def insertions(colour, fermions, tensor, prime):
    scalars = on_traceless(colour, tensor, prime)
    transposes = [x.T for x in scalars]
    left = [tensor[2]]
    right = [tensor[3]]
    for _ in range(3):
        left.append(alternate(transposes, left[-1], prime))
        right.append(alternate(scalars, right[-1], prime))
    assert pairing(left[3], right[0], prime) == pairing(left[0], right[3], prime)
    rows = []
    if fermions is None:
        commutators = [bracket(scalars[(i+1)%3], scalars[(i+2)%3], prime) for i in range(3)]
        triple = int(np.trace(scalars[0]@commutators[0]))%prime
        rows.append(pairing(left[0], right[3], prime)*triple%prime)
        assert left[0].size*(prime-1)**2 < 2**63
        for before in range(4):
            pair = np.zeros((4,4), dtype=np.int64)
            for x, comm in zip(scalars, commutators):
                lo = np.stack([act(comm.T,left[before],axis,prime).ravel() for axis in range(4)])
                hi = np.stack([act(x,right[3-before],axis,prime).ravel() for axis in range(4)])
                pair += (lo@hi.T)%prime
            rows.extend((pair%prime).ravel().tolist())
        return np.array(rows, dtype=np.int64)[:,None]
    for odd in fermions:
        psi = on_traceless(odd, tensor, prime)
        field = sum(bracket(x,b,prime) for x,b in zip(scalars,psi))%prime
        rows.append([pairing(left[before], act(field,right[3-before],axis,prime),prime)
                     for before in range(4) for axis in range(4)])
    return np.array(rows, dtype=np.int64).T


def calculate(work, length, workers, cache=None):
    prior, joined = [json.loads((HERE.parent/"f4"/n).read_bytes())
                     for n in ("primary.json","joined_images.json")]
    prime = prior["prime"]
    record = next(row for row in prior["records"] if row["length"]==length)
    saved = next(row for row in joined["records"] if row["length"]==length)
    if cache:
        old = np.load(cache, allow_pickle=False)
        assert m.minor(old,prime) == saved["minor"]
    else:
        old,_,_ = inherited(record,saved,work,workers,native(work,prime,"jordan"),
                            {27:m.library(work,27,prime)},prime)
    tensor = tensors(prime)
    blocks = []
    letters = LETTERS[:3 if length==11 else 6]
    for recipe in record["frames"]:
        f = frame(letters,recipe["points"],recipe["k"],recipe["degree"],prime,recipe["seed"])
        for point in range(recipe["points"]):
            colour = [f.L[l][1][point,0] for l in LETTERS[:3]]
            odd = None if length==11 else [
                [SIGNS[i]*f.L[LETTERS[i]][1][point,j]%prime for i in (3,4,5)]
                for j in range(recipe["k"])]
            blocks.append(insertions(colour,odd,tensor,prime))
        print("length",length,"points",len(blocks),flush=True)
    new = np.concatenate(blocks,axis=1)
    assert old.shape[1]==new.shape[1]
    extended = m.minor(np.vstack((old,new)),prime)
    return dict(prime=prime,length=length,previous=saved["minor"],augmented=extended,
                inserted_values=new.tolist(),new_rows=new.shape[0])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--cache",type=Path)
    parser.add_argument("--length",type=int,choices=(10,11),required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    args.work.mkdir(exist_ok=True,parents=True)
    data=calculate(args.work,args.length,args.workers,args.cache)
    args.out.write_text(json.dumps(data,separators=(",",":"))+"\n")
    print(data["previous"]["rank"], data["augmented"]["rank"])
