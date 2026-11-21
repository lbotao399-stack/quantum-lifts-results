"""Resolve the (s,s,1^c) Taylor representation using single-flavour classical chains."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess


def weights(spin, tail, upper_flavour):
    if spin < 1 or tail < 0:
        raise ValueError("expected two equal positive rows followed by a nonnegative tail")
    if upper_flavour:
        return [(tail + 2 + k, 0, 0, spin - 1, spin - 1 - k) for k in range(spin)]
    return [(tail - k, 0, 0, spin + k, spin) for k in range(tail + 1)]


def saved(path, record):
    temp = path.with_suffix(".part")
    temp.write_text(json.dumps(record, indent=2) + "\n")
    temp.replace(path)


def compute(engine, group, spin, tail, lengths, indices, upper_flavour, out,
            *, threads, workers, prime, seed):
    out.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(engine.read_bytes()).hexdigest()
    sectors = weights(spin, tail, upper_flavour)
    if len(set(indices)) != len(indices) or any(i < 0 or i >= len(sectors) for i in indices):
        raise ValueError("invalid sector indices")
    def one(index):
        weight = sectors[index]
        filename = f"{group}_{spin}_{tail}_{index}"
        campaign_seed = seed + 1009 * index
        args = [str(engine), "chain", group, *map(str, weight), *map(str, lengths),
                str(campaign_seed), str(prime)]
        environment = dict(os.environ, A_NORHO="1", A_VERB="1",
                           PAR_THREADS=str(max(1, threads // workers)))
        with (out.parent / (filename + ".log")).open("w") as log:
            raw = subprocess.check_output(args, stderr=log, env=environment)
        data = json.loads(raw)
        assert data["gauge"] == group and data["sector"] == list(weight)
        assert data["prime"] == prime and data["seed"] == campaign_seed
        assert [b["Y"] for b in data["blocks"]] == list(range(lengths[0], lengths[1] + 1))
        for b in data["blocks"]:
            assert b["dim"] - b["rank_in"] - b["rank_out"] == b["h"] >= 0
        for left, right in zip(data["blocks"], data["blocks"][1:]):
            assert left["rank_out"] == right["rank_in"]
        data.update(row_length=spin, tail=tail, upper_flavour=upper_flavour, sector_index=index,
                    engine_sha256=digest, arguments=args[1:])
        path = out / (filename + ".json")
        saved(path, data)
        print(group, weight, [(b["Y"], b["h"]) for b in data["blocks"]], flush=True)
        return path
    with ThreadPoolExecutor(workers) as pool:
        return list(pool.map(one, indices))


def resolve(records):
    first = records[0]
    group, spin, tail, direction, prime = (
        first[k] for k in ("gauge", "row_length", "tail", "upper_flavour", "prime"))
    sectors = weights(spin, tail, direction)
    ordered = sorted(records, key=lambda r: r["sector_index"])
    assert [r["sector_index"] for r in ordered] == list(range(len(sectors)))
    for weight, row in zip(sectors, ordered, strict=True):
        assert (row["sector"], row["gauge"], row["row_length"], row["tail"], row["upper_flavour"],
                row["prime"]) == (list(weight), group, spin, tail, direction, prime)
    blocks = []
    for terms in zip(*(r["blocks"] for r in ordered), strict=True):
        y = terms[0]["Y"]
        assert all(b["Y"] == y for b in terms)
        value = {key: sum((-1)**i * b[key] for i, b in enumerate(terms))
                 for key in ("dim", "rank_in", "rank_out", "h")}
        assert min(value.values()) >= 0
        assert value["dim"] - value["rank_in"] - value["rank_out"] == value["h"]
        blocks.append(dict(Y=y, **value))
    return dict(group=group, partition=[spin, spin] + [1] * tail, prime=prime,
                upper_flavour=direction, blocks=blocks)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--group", required=True)
    p.add_argument("--spin", type=int, required=True, help="common length of the first two rows")
    p.add_argument("--tail", type=int, required=True)
    p.add_argument("--lengths", type=int, nargs=2, required=True)
    p.add_argument("--indices", type=int, nargs="+")
    p.add_argument("--lower-flavour", action="store_true")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--prime", type=int, default=2147483647)
    p.add_argument("--seed", type=int, default=876331)
    a = p.parse_args()
    direction = not a.lower_flavour
    indices = a.indices if a.indices is not None else list(range(len(weights(a.spin, a.tail, direction))))
    files = compute(a.engine.resolve(), a.group, a.spin, a.tail, a.lengths, indices,
                    direction, a.out, threads=a.threads, workers=a.workers, prime=a.prime, seed=a.seed)
    if len(files) == len(weights(a.spin, a.tail, direction)):
        saved(a.out / "resolved.json", resolve([json.loads(f.read_bytes()) for f in files]))
