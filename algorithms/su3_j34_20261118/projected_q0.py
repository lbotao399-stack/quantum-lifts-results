"""Compute joint minors of classical images in selected level-34 target profiles."""
import argparse
import gzip
import hashlib
import json
import math
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import chain
from pathlib import Path

import numpy as np
from flint import nmod_mat

HERE = Path(__file__).resolve().parent
RUNTIME = HERE.parent / "su3_independent_prime_20261118/runtime/work/sp3-native.nosync"


def pattern(row):
    return tuple((tuple(key), count) for key, count in row["profile"])


def setup(runtime, plan, targets, sources, prime, narrow=False):
    global algebra, choices, charge, evaluate, frames, project, random_word, source_rows, target_rows, modulus
    sys.path.insert(0, str(runtime))
    import su_trace_adapter_v1 as algebra
    from g2_full_classical_screen_v1 import contents, sector
    from fortuity.full.evalg import Frame
    from native_su_projection_v1 import project
    from native_su_profile_derivation_v1 import evaluate
    from sp_full_profiles import profile
    from su3_stream_rank_v3 import random_word
    weight = tuple(plan["mu"])
    length = plan["source_length"]
    assert plan["parity"] == 1
    charge = sector(weight, length)
    source_rows = [plan["sources"][i] for i in sources]
    target_rows = [plan["targets"][i] for i in targets]
    modulus = prime
    contents_by_profile = {}
    for content in contents(weight, length):
        contents_by_profile.setdefault(profile(content), []).append(content)
    choices = [contents_by_profile[pattern(s)] for s in source_rows]
    # The native differential needs the union of the predecessor and image alphabets.
    alphabet = sorted({letter for size in (length, length+1) for c in contents(weight, size) for letter in c})
    frames = []
    for t in target_rows:
        degree = t["fermions"]
        generators = degree+2
        frame_dimension = min(sum(s["dimension"] for s in source_rows), t["dimension"]) if narrow else t["dimension"]
        points = math.ceil((frame_dimension+8)/math.comb(generators, degree))
        frames.append(Frame(3, points, generators, [degree], alphabet,
                            seed=2026111801+7919*t["index"], p=prime))


def decode(word):
    return tuple(tuple((a,b,tuple(theta)) for a,b,theta in trace) for trace in word)


def values(selected):
    result = np.zeros((len(selected), sum(f.ncols for f in frames)), dtype=np.int64)
    for i in sorted({i for i,_ in selected}):
        positions = [j for j,(s,_) in enumerate(selected) if s == i]
        expr = [project(algebra.parity_project({selected[j][1]:1},1), charge, modulus) for j in positions]
        start = 0
        for t,frame in zip(target_rows,frames):
            stop = start+frame.ncols
            if t["index"] in source_rows[i]["targets"]:
                result[positions,start:stop] = evaluate(expr, frame, pattern(source_rows[i]), pattern(t))
            start = stop
    return result


def sample(batch):
    rng = random.Random(2026111809+104729*batch)
    selected = []
    for j in range(8):
        i = (8*batch+j) % len(choices)
        for _ in range(128):
            word = random_word(rng.choice(choices[i]),rng,6)
            if word is not None:
                selected.append((i,word))
                break
    return selected,values(selected)


def write_json(path, data):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data,indent=2)+"\n")
    tmp.replace(path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runtime", type=Path, default=RUNTIME)
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--direction", choices=("incoming","outgoing"), default="incoming")
    p.add_argument("--targets", type=int, nargs="+", required=True)
    p.add_argument("--sources", type=int, nargs="+")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--prime", type=int, default=1000003,choices=(1000003,1000033))
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--seconds", type=int, default=600)
    p.add_argument("--batches", type=int, default=20000)
    p.add_argument("--replay", type=Path)
    p.add_argument("--narrow", action="store_true", help="Size a single target frame by the smaller cochain dimension.")
    args = p.parse_args()
    record = json.loads(args.plan.read_text())
    plan = record[args.direction]
    assert len(set(args.targets)) == len(args.targets)
    possible = [s["index"] for s in plan["sources"] if set(args.targets).intersection(s["targets"])]
    sources = possible if args.sources is None else args.sources
    assert len(set(sources)) == len(sources) and set(sources) <= set(possible)
    assert not args.narrow or len(args.targets) == 1
    init = (args.runtime.resolve(),plan,args.targets,sources,args.prime,args.narrow)
    setup(*init)
    from su_modular_rowspace_v1 import RowSpace
    width = sum(f.ncols for f in frames)
    source_dimension = sum(s["dimension"] for s in source_rows)
    target_dimension = sum(t["dimension"] for t in target_rows)
    identity = dict(order=record["order"], prime=args.prime, sources=sources, targets=args.targets,
                    plan_sha256=hashlib.sha256(args.plan.read_bytes()).hexdigest(), frame_width=width,
                    source_dimension=source_dimension,target_dimension=target_dimension)
    if args.narrow:
        identity["frame_rule"] = "min(source_dimension,target_dimension)+8 before rounding to points"
    if args.direction == "outgoing":
        identity["direction"] = args.direction
    if args.replay:
        cert = json.loads(gzip.decompress(args.replay.read_bytes()))
        assert cert["map"] == identity
        recipes = [(i,decode(w)) for i,w in cert["rows"]]
        assert len(set(recipes)) == len(set(cert["columns"])) == len(recipes)
        blocks = [values(recipes[j:j+8])[:,cert["columns"]] for j in range(0,len(recipes),8)]
        det = int(nmod_mat(np.vstack(blocks).tolist(),args.prime).det())
        assert det == cert["determinant"] and det
        print(json.dumps(dict(order=len(recipes),determinant=det,prime=args.prime)),flush=True)
        return
    args.output.mkdir(parents=True,exist_ok=True)
    space = RowSpace(width,args.prime)
    rows, recipes = [], []
    batch = 0
    clock = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers,initializer=setup,initargs=init) as pool:
        while (batch < args.batches and time.monotonic()-clock < args.seconds
               and len(recipes) < min(source_dimension,target_dimension)):
            stop = min(args.batches,batch+4*args.workers)
            packets = list(pool.map(sample,range(batch,stop)))
            words = list(chain.from_iterable(w for w,_ in packets))
            block = np.vstack([v for _,v in packets])
            keep = space.extend(block)
            recipes.extend(words[j] for j in keep)
            rows.extend(block[j].astype(np.uint32) for j in keep)
            assert len(recipes) <= min(source_dimension,target_dimension)
            batch = stop
            print(json.dumps(dict(batches=batch,rank=len(recipes))),flush=True)
    cols = list(map(int,space.pivots))
    if not rows:
        raise RuntimeError("No nonzero pivot was obtained")
    minor = np.asarray(rows,dtype=np.uint32)[:,cols]
    det = int(nmod_mat(minor.tolist(),args.prime).det())
    assert det and len(set(recipes)) == len(cols) == len(recipes)
    cert = dict(map=identity,rows=recipes,columns=cols,determinant=det)
    payload = gzip.compress(json.dumps(cert,separators=(",",":")).encode(),mtime=0)
    target = args.output/"minor.json.gz"
    tmp = target.with_suffix(".tmp")
    tmp.write_bytes(payload)
    tmp.replace(target)
    result = dict(**identity,batches=batch,rank_lower_bound=len(cols),rank_ceiling=min(source_dimension,target_dimension),
                  determinant=det,sha256=hashlib.sha256(payload).hexdigest(),bytes=len(payload))
    write_json(args.output/"result.json",result)
    print(json.dumps(result),flush=True)


if __name__ == "__main__":
    main()
