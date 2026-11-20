"""Classical singlet cohomology before rank-three trace identities, grouped by trace degrees."""
import argparse
import gzip
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mixed_incoming as m
from mixed33_exact import apply, exact, sha
from mixed_middle import build


def degree(w):
    return tuple(sorted((sum(a+b for a,b,theta in tr), sum(len(theta) for a,b,theta in tr))
                        for tr in w))


def decode_rows(data, vector_name, dictionary_name):
    words = list(map(m.decode, data[dictionary_name]))
    return [{words[i]: c for i,c in row} for row in data[vector_name]]


def read_plan(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def compute(task):
    label, incoming, currents, outgoing, path = task
    if path.is_file():
        return
    domain = sorted({w for row in incoming+currents for w in row})
    codomain = sorted({w for row in outgoing for w in row})
    def packed(rows, words):
        ids = {w:i for i,w in enumerate(words)}
        return [{ids[w]:v for w,v in row.items()} for row in rows]
    inc = packed(incoming, domain)
    src = packed(currents, domain)
    out = packed(outgoing, codomain)
    independent, _, _ = exact(inc)
    active, _, kernel = exact(out)
    assert all(not apply(v,out) for v in kernel)
    cycles = [apply(v,src) for v in kernel]
    joint, _, _ = exact([inc[i] for i in independent] + cycles)
    keep = [i-len(independent) for i in joint if i >= len(independent)]
    assert len(keep) == len(src)-len(active)-len(independent)
    results = [[[m.encode(domain[j]),c] for j,c in sorted(cycles[i].items())] for i in keep]
    m.write(path, dict(trace_degrees=label, source_singlets=len(src), incoming_rank=len(independent),
                       outgoing_rank=len(active), representatives=results))


def run(t, work, first_plan, workers):
    work.mkdir(parents=True, exist_ok=True)
    second = build(t, work, workers, length=11-t)
    ip,op = read_plan(first_plan),read_plan(second)
    assert ip["mu"] == op["mu"] and ip["length"]+1 == op["length"]
    inc, src, out = (decode_rows(doc,k,d) for doc,k,d in
                     ((ip,"images","target_dictionary"),(op,"sources","source_dictionary"),
                      (op,"images","target_dictionary")))
    grouped = defaultdict(lambda: ([],[],[]))
    for row in inc:
        if row:
            labs = {degree(w) for w in row}
            assert len(labs) == 1
            grouped[labs.pop()][0].append(row)
    for x,y in zip(src,out):
        labs = {degree(w) for w in x.keys() | y.keys()}
        assert len(labs) == 1
        _,sources,images = grouped[labs.pop()]
        sources.append(x)
        images.append(y)
    tasks = [(k,*parts,work/f"complex-{i:04}.json.gz") for i,(k,parts) in enumerate(sorted(grouped.items()))]
    print("formal blocks",len(tasks),"largest", max(len(t[1])+len(t[2]) for t in tasks),flush=True)
    with ProcessPoolExecutor(workers) as pool:
        for _ in pool.map(compute, sorted(tasks,key=lambda v:-len(v[1])-len(v[2]))):
            pass
    records=[read_plan(row[-1]) for row in tasks]
    n=sum(len(row["representatives"]) for row in records)
    m.write(work/"classical.json.gz", dict(mu=op["mu"], length=op["length"],
            incoming_plan_sha256=sha(first_plan), outgoing_plan_sha256=sha(second),
            source_singlets=len(src), formal_classical_dimension=n, blocks=records))
    print("formal cohomology",n,"word counts",[len(r) for b in records for r in b["representatives"]],flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t",type=int,choices=(1,2,3),required=True)
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--first-plan",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=8)
    assert __debug__
    run(**vars(parser.parse_args()))
