"""Checkpoint disjoint ranges of a complete scalar-chart expansion and verify coverage."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path

from su3_subregular_image import expand_words, source, write


def digest(data):
    return hashlib.sha256(json.dumps(data,separators=(",",":")).encode()).hexdigest()


def part(task):
    start, words, degree, folder = task
    path = folder/f"piece-{start:07}.json.gz"
    metadata = dict(start=start, stop=start+len(words), degree=degree, input_sha256=digest(words))
    if path.is_file():
        stored = json.loads(gzip.decompress(path.read_bytes()))
        if stored["metadata"] != metadata or digest(stored["coefficients"]) != stored["output_sha256"]:
            raise ValueError(f"incorrect checkpoint {path}")
    else:
        values = expand_words(words, degree)
        terms = [[w.hex(),n] for w,n in sorted(values.items())]
        tmp = path.with_suffix(".partial.gz")
        write(tmp,dict(metadata=metadata,output_sha256=digest(terms),coefficients=terms))
        tmp.replace(path)
    return path


def build(level, degree, source_path, work, workers, chunk):
    work.mkdir(exist_ok=True,parents=True)
    packet = source(level,source_path,degree)
    selected = packet.pop("selected")
    tasks = ((i,selected[i:i+chunk],degree,work) for i in range(0,len(selected),chunk))
    with ProcessPoolExecutor(workers) as pool:
        paths = list(pool.map(part,tasks))
    total = Counter()
    stop = 0
    for path in paths:
        stored=json.loads(gzip.decompress(path.read_bytes()))
        metadata=stored["metadata"]
        assert metadata["start"]==stop and metadata["input_sha256"]==digest(selected[stop:metadata["stop"]])
        assert digest(stored["coefficients"])==stored["output_sha256"]
        stop=metadata["stop"]
        total.update({bytes.fromhex(w):n for w,n in stored["coefficients"]})
    assert stop==len(selected)
    degrees={}
    for k in range(degree+1):
        terms=[[w.hex(),n] for w,n in sorted((w,n) for w,n in total.items() if len(w)==k and n)]
        if terms:
            path=work/f"degree-{k}.json.gz"
            write(path,terms)
            degrees[str(k)]=dict(terms=len(terms),payload_sha256=digest(terms),file=path.name)
    receipt=dict(**packet, selected_words=len(selected), parts=len(paths), homogeneous=degrees)
    write(work/"verified.json",receipt)
    return receipt


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level",type=int,choices=(32,34),required=True)
    p.add_argument("--degree",type=int,required=True)
    p.add_argument("--input",type=Path)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--chunk",type=int,default=1800)
    args=p.parse_args()
    if not __debug__ or min(args.workers,args.chunk)<1 or args.degree<0:
        p.error("positive worker and chunk counts, nonnegative degree, and assertions are required")
    print(build(args.level,args.degree,args.input,args.work,args.workers,args.chunk),flush=True)
