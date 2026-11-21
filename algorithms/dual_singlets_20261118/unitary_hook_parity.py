"""Charge conjugation on the classical (s,s,1^c) Taylor component of SU(N)."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from taylor_hook_cohomology import saved, weights
from trace_descent import dimensions


def calculate(engine, colours, spin, tail, length, prime, seed, indices, directory, threads):
    directory.mkdir(parents=True, exist_ok=True)
    environment = dict(os.environ, PAR_THREADS=str(threads))
    sectors = weights(spin, tail, True)
    if not indices or len(set(indices)) != len(indices) or any(i not in range(spin) for i in indices):
        raise ValueError("invalid horizontal-strip indices")
    for index in indices:
        weight = sectors[index]
        args = [str(engine.resolve()), f"SU{colours}", *map(str, weight),
                str(length), str(seed + 1009*index), str(prime)]
        record = json.loads(subprocess.check_output(args, env=environment))
        checks = (record["gauge"], record["sector"], record["prime"], record["length"], record["seed"])
        if checks != (f"SU{colours}", list(weight), prime, length, seed+1009*index):
            raise ValueError("engine record does not match requested sector")
        character_weight = (*weight[3:], *weight[:3])
        for parity, group in [("even", "SU+"), ("odd", "SU-")]:
            dimension = sum(dimensions(character_weight, length, colours, group).values())
            record[parity+"_dim"] = dimension
            h = dimension - record[parity+"_in"] - record[parity+"_out"]
            if h < 0:
                raise ValueError("adjacent image ranks exceed character multiplicity")
            record[parity+"_h"] = h
        if record["even_dim"] + record["odd_dim"] != record["cochains"]:
            raise ValueError("independent colour characters disagree")
        record.update(row_length=spin, tail=tail, sector_index=index)
        saved(directory / f"sector_{index}.json", record)
        print("Sector", index, "classical dimensions", record["even_h"], record["odd_h"], flush=True)


def resolve(records):
    first = records[0]
    spin, tail = first["row_length"], first["tail"]
    prime, group, length = first["prime"], first["gauge"], first["length"]
    ordering = sorted(records, key=lambda row: row["sector_index"])
    if [r["sector_index"] for r in ordering] != list(range(spin)):
        raise ValueError("incomplete or duplicate horizontal strips")
    for weight, row in zip(weights(spin, tail, True), ordering, strict=True):
        if (row["sector"], row["gauge"], row["prime"], row["length"], row["row_length"], row["tail"]) != (
                list(weight), group, prime, length, spin, tail):
            raise ValueError("incompatible Taylor sectors")
    blocks = {}
    for parity in ("even", "odd"):
        values = {key:sum((-1)**i*r[parity+"_"+key] for i, r in enumerate(ordering))
                  for key in ("dim", "in", "out", "h")}
        if min(values.values()) < 0 or values["h"] != values["dim"]-values["in"]-values["out"]:
            raise ValueError("inconsistent Taylor multiplicities")
        blocks[parity] = values
    return dict(group=group, partition=[spin, spin]+[1]*tail, length=length, prime=prime, parities=blocks)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--colours", type=int, choices=(2, 3), default=3)
    p.add_argument("--spin", type=int, default=5)
    p.add_argument("--tail", type=int, default=3)
    p.add_argument("--length", type=int, default=9)
    p.add_argument("--prime", choices=(1000003,1000033,2147483647), type=int, default=2147483647)
    p.add_argument("--seed", type=int, default=375901)
    p.add_argument("--indices", type=int, nargs="+")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args()
    indices = a.indices or list(range(a.spin))
    calculate(a.engine, a.colours, a.spin, a.tail, a.length, a.prime, a.seed, indices, a.out, a.threads)
    if sorted(indices) == list(range(a.spin)):
        saved(a.out/"resolved.json", resolve([json.loads((a.out/f"sector_{i}.json").read_bytes()) for i in indices]))
