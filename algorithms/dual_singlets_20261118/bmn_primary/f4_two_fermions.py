"""Apply Q0 to two scalar--fermion insertions in the four-arm Albert contraction."""
import argparse
import json
from pathlib import Path
import numpy as np

from f4_cubic_boundaries import HERE, LETTERS, SIGNS, bracket, frame, m
from f4_scalar_nine import act, alternate, on_traceless, tensors


SLOTS = ((0,0),(1,1),(0,1),(1,0),(0,3),(3,0))


def paired(value, first, second, transpose, prime):
    out = np.zeros((len(SLOTS),value.size),dtype=np.int64)
    for x,y in zip(first,second):
        before = y.T if transpose else x
        after = x.T if transpose else y
        done = {}
        for row,(a,b) in enumerate(SLOTS):
            if transpose:
                a,b=b,a
            if a not in done:
                done[a]=act(before,value,a,prime)
            out[row]+=act(after,done[a],b,prime).ravel()
    return out%prime


def point_rows(scalars, odd, tensor, prime):
    right=[tensor[3]]
    left=[tensor[2]]
    for _ in range(2):
        right.append(alternate(scalars,right[-1],prime))
        left.append(alternate([x.T for x in scalars],left[-1],prime))
    comm=[bracket(scalars[(i+1)%3],scalars[(i+2)%3],prime) for i in range(3)]
    blocks=[]
    assert left[0].size*(prime-1)**2 < 2**63
    for before in range(3):
        closed_left=paired(left[before],scalars,comm,True,prime)
        closed_right=paired(right[2-before],scalars,comm,False,prime)
        columns=[]
        for psi in odd:
            odd_left=paired(left[before],scalars,psi,True,prime)
            odd_right=paired(right[2-before],scalars,psi,False,prime)
            values=(closed_left@odd_right.T-odd_left@closed_right.T)%prime
            columns.append(values.ravel())
        blocks.append(np.stack(columns,axis=1))
    return np.vstack(blocks)


def calculate(cache,columns):
    from f4_bmn_double_jordan import inherited
    from f4_joined_boundaries import native
    from tempfile import TemporaryDirectory
    prior,joined=(json.loads((HERE.parent/"f4"/n).read_bytes())
                  for n in ("primary.json","joined_images.json"))
    rec,base=(next(r for r in f["records"] if r["length"]==10) for f in (prior,joined))
    prime=prior["prime"]
    if cache:
        old=np.load(cache,allow_pickle=False)
    else:
        with TemporaryDirectory() as directory:
            work=Path(directory)
            old,_,_=inherited(rec,base,work,4,native(work,prime,"jordan"),
                              {27:m.library(work,27,prime)},prime)
    assert m.minor(old,prime)==base["minor"]
    if not 80<=columns<=old.shape[1]:
        raise ValueError("The number of mixed evaluation columns is outside the saved frame")
    old=old[:,:columns]
    tensor=tensors(prime)
    blocks=[]
    for row in rec["frames"]:
        f=frame(LETTERS[:6],row["points"],row["k"],1,prime,row["seed"])
        for pt in range(row["points"]):
            scalars=on_traceless([f.L[l][1][pt,0] for l in LETTERS[:3]],tensor,prime)
            odd=[on_traceless([SIGNS[i]*f.L[LETTERS[i]][1][pt,j]%prime for i in (3,4,5)],
                              tensor,prime) for j in range(row["k"])]
            blocks.append(point_rows(scalars,odd,tensor,prime))
            if len(blocks)*row["k"]>=columns:
                break
        print("points",len(blocks),flush=True)
        if len(blocks)*row["k"]>=columns:
            break
    new=np.concatenate(blocks,axis=1)[:,:columns]
    return dict(prime=prime,slots=SLOTS,columns=columns,previous=m.minor(old,prime),
                augmented=m.minor(np.vstack((old,new)),prime),inserted_values=new.tolist())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path)
    p.add_argument("--columns",type=int,default=81)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=calculate(a.cache,a.columns)
    a.out.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(result["previous"]["rank"],result["augmented"]["rank"])
