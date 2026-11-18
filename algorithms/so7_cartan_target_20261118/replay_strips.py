"""Parallel row reconstruction and checkpointed determinant replay for the SO(7) minor."""
import argparse
import gzip
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from flint import nmod_mat

from incoming_minor import evaluate_batch
from incoming_rows import setup
from reduce_incoming import P, product


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def build(cert, args):
    panels=cert["panels"]
    columns=[j for b in panels for j in b["columns"]]
    rows=[row for b in panels for row in b["words"]]
    order=sorted(range(len(rows)),key=lambda j:rows[j][0])
    init=tuple(str(p.resolve()) for p in (args.runtime,args.native,args.plan,args.functional))
    folder=args.output/"strips"
    folder.mkdir(parents=True,exist_ok=True)
    expected=0
    with ProcessPoolExecutor(max_workers=args.workers,initializer=setup,initargs=init) as pool:
        for start in range(0,len(rows),32*args.shards):
            pos=start+32*args.shard
            idx=order[pos:pos+32]
            if not idx:
                continue
            expected+=len(idx)
            path=folder/f"rows_{pos:05}.npz"
            if path.exists():
                with np.load(path,allow_pickle=False) as z:
                    assert z["matrix"].shape==(len(idx),len(columns)) and np.array_equal(z["indices"],idx)
                continue
            tasks=[([rows[j] for j in idx[k:k+2]],cert["widths"],columns) for k in range(0,len(idx),2)]
            block=np.vstack(list(pool.map(evaluate_batch,tasks))).astype("<u4")
            assert block.shape==(len(idx),len(columns)) and int(block.max())<P
            temp=path.with_suffix(".tmp")
            with temp.open("wb") as f:
                np.savez(f,matrix=block,indices=np.asarray(idx,dtype="<u4"))
            temp.replace(path)
            print(json.dumps(dict(shard=args.shard,rows=expected)),flush=True)
    print(json.dumps(dict(shard=args.shard,complete=True,rows=expected)),flush=True)


def verify(cert, args):
    panels=cert["panels"]
    n=sum(len(p["columns"]) for p in panels)
    matrix=np.empty((n,n),dtype="<u4")
    used=set()
    for path in sorted((args.output/"strips").glob("rows_*.npz")):
        with np.load(path,allow_pickle=False) as z:
            indices=list(map(int,z["indices"]))
            assert indices and min(indices)>=0 and max(indices)<n
            assert len(set(indices))==len(indices) and used.isdisjoint(indices)
            used.update(indices)
            matrix[indices]=z["matrix"]
    assert len(used)==n==cert["rank"]
    earlier=[]
    start=0
    det=1
    for i,panel in enumerate(panels):
        stop=start+len(panel["columns"])
        values=matrix[start:stop].astype(np.int64)
        for s,t,upper in earlier:
            values=(values-product(values[:,s:t],upper))%P
        square=nmod_mat(values[:,start:stop].tolist(),P)
        d=int(square.det())
        assert d==panel["determinant"] and d
        det=det*d%P
        inverse=np.asarray(square.inv().tolist(),dtype=np.int64)
        earlier.append((start,stop,product(inverse,values).astype(np.float64)))
        start=stop
        print(json.dumps(dict(panel=i,rank=stop)),flush=True)
    assert det==cert["determinant"] and start==n
    record=dict(prime=P,order=n,panels=len(panels),determinant=det,
                recomputed_entries=n*n,recipes_sha256=digest(args.archive),
                matrix_sha256=hashlib.sha256(memoryview(matrix)).hexdigest())
    (args.output/"incoming_replay.json").write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps(record),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--verify",action="store_true")
    p.add_argument("--runtime",type=Path)
    p.add_argument("--native",type=Path)
    p.add_argument("--plan",type=Path)
    p.add_argument("--functional",type=Path)
    p.add_argument("--shard",type=int,default=0)
    p.add_argument("--shards",type=int,default=1)
    p.add_argument("--workers",type=int,default=12)
    a=p.parse_args()
    assert digest(a.archive)=="6fbdb485472d3e9c6385e15a8d8bfc2cf690642a00f15dfb5e96321f44a4ed20"
    data=json.loads(gzip.decompress(a.archive.read_bytes()))
    assert data["prime"]==P and 0<=a.shard<a.shards
    (verify if a.verify else build)(data,a)
