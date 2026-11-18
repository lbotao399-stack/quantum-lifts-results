"""Extract and recompute a square minor of the simultaneous closure and pairing matrix."""
import argparse
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from itertools import islice
from pathlib import Path

import numpy as np
from flint import nmod_mat
from reduce_incoming import Echelon, P, product
import incoming_rows

HERE = Path(__file__).resolve().parent
PROFILE_ORDER = (2,0,4,5,1,9,8,7,3,6)


def offsets(widths):
    starts = []
    for w in widths:
        starts.append(sum(widths[:len(starts)]))
    return starts, sum(widths)


def joined(values, widths):
    starts, width = offsets(widths)
    output = np.zeros((len(values["ell_Q1"]),width+1),dtype=np.uint32)
    for key, data in values.items():
        if key != "ell_Q1":
            j=int(key)
            assert data.shape[1]==widths[j]
            output[:,starts[j]:starts[j]+widths[j]]=data
    output[:,-1]=values["ell_Q1"]
    return output


def determinant(panels):
    value=1
    for block in panels:
        value=value*block["determinant"]%P
    return value


def extract(packed, raw, output):
    coverage=json.loads((packed/"coverage.json").read_text())
    widths=[next(r["widths"][str(i)] for r in coverage if str(i) in r["widths"]) for i in range(3)]
    last=sum(widths)
    space=Echelon(last+1,minors=True)
    tags={}
    for i in PROFILE_ORDER:
        record=coverage[i]
        assert record["index"]==i
        tags[i]=np.load(packed/("p%02d_tags.npy"%i),mmap_mode="r")
        arrays={k:np.load(packed/("p%02d_%s.npy"%(i,k)),mmap_mode="r")
                for k in record["widths"] if k!="source"}
        arrays["ell_Q1"]=np.load(packed/("p%02d_rhs.npy"%i),mmap_mode="r")
        for j in range(0,record["rows"],128):
            stop=min(j+128,record["rows"])
            block=joined({k:v[j:stop] for k,v in arrays.items()},widths)
            space.add(block,[(i,s) for s in range(j,stop)])
            if space.rank>=12526:
                break
        print(json.dumps(dict(profile=i,rank=space.rank)),flush=True)
        if space.rank>=12526:
            break
    assert space.rank==12526
    panels=[]
    for panel in space.panels:
        words=[]
        for i,j in panel["ids"]:
            serial,sub=map(int,tags[i][j])
            path=raw/("p%02d_b%07d.json.gz"%(i,serial))
            words.append([i,json.loads(gzip.decompress(path.read_bytes()))[sub]])
        panels.append(dict(columns=list(map(int,panel["pivots"])),words=words,determinant=panel["determinant"]))
    columns=[j for panel in panels for j in panel["columns"]]
    assert len(set(columns))==len(columns)==space.rank and last in columns
    cert=dict(prime=P,widths=widths,rank=space.rank,determinant=determinant(panels),panels=panels)
    assert cert["determinant"]
    payload=gzip.compress(json.dumps(cert,separators=(",",":")).encode(),mtime=0)
    output.mkdir(parents=True,exist_ok=True)
    (output/"incoming_minor.json.gz").write_bytes(payload)
    info=dict(prime=P,order=space.rank,determinant=cert["determinant"],
              panels=len(panels),pairing_column=last,bytes=len(payload),
              recipes_sha256=hashlib.sha256(payload).hexdigest())
    (output/"minor_result.json").write_text(json.dumps(info,indent=2)+"\n")
    print(json.dumps(info),flush=True)


def evaluate_batch(task):
    rows,widths,columns=task
    out=np.zeros((len(rows),len(columns)),dtype=np.uint32)
    for i in sorted({i for i,_ in rows}):
        positions=[j for j,(idx,_) in enumerate(rows) if idx==i]
        words=[tuple(tuple((u,v,tuple(s)) for u,v,s in t) for t in rows[j][1]) for j in positions]
        values=incoming_rows.evaluate(i,words,False)
        out[positions]=joined(values,widths)[:,columns]
    return out


def replay(archive, runtime, native, plan, functional, output, workers):
    cert=json.loads(gzip.decompress(archive.read_bytes()))
    widths, blocks=cert["widths"],cert["panels"]
    columns=[j for panel in blocks for j in panel["columns"]]
    assert cert["prime"]==P and len(columns)==len(set(columns))==cert["rank"] and sum(widths) in columns
    offset=0
    earlier=[]
    init=tuple(str(p.resolve()) for p in (runtime,native,plan,functional))
    with ProcessPoolExecutor(max_workers=workers,initializer=incoming_rows.setup,initargs=init) as pool:
        for number,panel in enumerate(blocks):
            count=len(panel["words"])
            assert count==len(panel["columns"])
            it=iter(panel["words"])
            groups=iter(lambda:list(islice(it,4)),[])
            values=np.vstack(list(pool.map(evaluate_batch,((group,widths,columns) for group in groups)))).astype(np.int64)
            for s,t,matrix in earlier:
                values=(values-product(values[:,s:t],matrix))%P
            stop=offset+count
            square=nmod_mat(values[:,offset:stop].tolist(),P)
            assert int(square.det())==panel["determinant"] and panel["determinant"]
            inverse=np.asarray(square.inv().tolist(),dtype=np.int64)
            earlier.append((offset,stop,product(inverse,values).astype(np.float64)))
            offset=stop
            print(json.dumps(dict(panel=number,verified_rank=offset)),flush=True)
    result=dict(prime=P,order=offset,determinant=determinant(blocks),
                recipes_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    assert result["determinant"]==cert["determinant"] and offset==cert["rank"]
    output.mkdir(parents=True,exist_ok=True)
    (output/"incoming_replay.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--packed",type=Path)
    p.add_argument("--raw",type=Path)
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--replay",type=Path)
    p.add_argument("--runtime",type=Path)
    p.add_argument("--native",type=Path)
    p.add_argument("--plan",type=Path,default=HERE/"source_profiles.json")
    p.add_argument("--functional",type=Path)
    p.add_argument("--workers",type=int,default=1)
    a=p.parse_args()
    if a.replay:
        replay(a.replay,a.runtime,a.native,a.plan,a.functional,a.output,a.workers)
    else:
        extract(a.packed,a.raw,a.output)
