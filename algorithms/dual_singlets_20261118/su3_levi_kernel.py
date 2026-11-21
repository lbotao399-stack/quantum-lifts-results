"""Compute and replay classical Cartan-kernel bounds on the full rank-one Levi of SU(3)."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import json
import os
from pathlib import Path
import subprocess

from su3_levi_branching import CASES, plan


def key(weight, length):
    return "_".join(map(str, (*weight, length)))


def seed_for(weight, length):
    seed = 104729
    for degree in (*weight, length):
        seed = 17 * seed + degree
    return seed


def check(data, weight, length, prime, gauge=None):
    mode = "classical_Cartan_kernel" if gauge else "SU2_classical_Cartan_kernel"
    if (data["gauge"], data["mode"], data["sector"], data["length"], data["prime"], data["seed"]) != (
            gauge or "SU2", mode, list(weight), length, prime, seed_for(weight, length)):
        raise ValueError("sector, prime or seed disagrees")
    size = data["cochains"]
    if not size:
        return 0, 0, 0
    sizes = {"incoming": (size, data["parent"]), "outgoing": (data["next"], size),
             "outgoing_and_Cartan": (data["stacked_rows"], size)}
    for name, (height, width) in sizes.items():
        minor = data[name]
        rank = minor["rank"]
        if not (0 <= rank <= min(height, width) and 0 < minor["det"] < prime):
            raise ValueError(f"invalid {name} minor at {weight}, {length}")
        for label, bound in (("rows", height), ("columns", width)):
            indices = minor[label]
            if (len(indices) != rank or len(set(indices)) != rank or
                    any(not 0 <= index < bound for index in indices)):
                raise ValueError(f"invalid {name} {label}")
    h = size - data["incoming"]["rank"] - data["outgoing"]["rank"]
    invisible = size - data["incoming"]["rank"] - data["outgoing_and_Cartan"]["rank"]
    if not 0 <= invisible <= h:
        raise ValueError("incompatible classical ranks")
    return size, h, invisible


def gzip_json(path, data):
    raw = (json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n").encode()
    tmp = path.with_suffix(".part")
    tmp.write_bytes(gzip.compress(raw, mtime=0))
    tmp.replace(path)


def compute(cases, engine, folder, workers, prime, gauge=None):
    folder.mkdir(parents=True, exist_ok=True)
    jobs = {(tuple(row["weight"]), row["length"]) for case in cases.values() for row in case["active"]
            if row["length"]}
    def run(task):
        weight, length = task
        path = folder / (key(weight, length) + ".json")
        if path.is_file():
            record = json.loads(path.read_bytes())
        else:
            args = [str(engine.resolve()), *map(str, weight), str(length), str(prime), str(seed_for(weight, length))]
            if gauge:
                args.append(gauge)
            record = json.loads(subprocess.check_output(args, env=dict(os.environ, PAR_THREADS="1")))
        check(record, weight, length, prime, gauge)
        if not path.is_file():
            tmp = path.with_suffix(".part")
            tmp.write_text(json.dumps(record, separators=(",", ":")) + "\n")
            tmp.replace(path)
        return key(weight, length), record

    with ThreadPoolExecutor(workers) as pool:
        tasks = [pool.submit(run, t) for t in sorted(jobs, key=lambda x: (-sum(x[0]), -x[1], x[0]))]
        out = {}
        for future in as_completed(tasks):
            label, result = future.result()
            out[label] = result
            if len(out) % 200 == 0:
                print(len(out), "/", len(jobs), flush=True)
    return out


def resolve(cases, records, prime):
    from trace_descent import dimensions
    output = {}
    for name, case in cases.items():
        totals = [0, 0, 0]
        remaining = []
        for row in case["active"]:
            weight, length, copies = row["weight"], row["length"], row["copies"]
            if length:
                counts = check(records[key(weight, length)], weight, length, prime)
            else:
                if any(weight):
                    raise ValueError("nontrivial weight assigned to the unit")
                counts = (1, 1, 0)
            for i, n in enumerate(counts):
                totals[i] += copies * n
            if counts[2]:
                remaining.append(dict(active_weight=weight, active_length=length, copies=copies,
                                      kernel_upper_per_copy=counts[2]))
        weight = case["weight"]
        colour = "U+" if case["colour_parity"] == 1 else "U-"
        direct = sum(dimensions((*weight[3:], *weight[:3]), case["length"], 2, colour).values())
        if direct != totals[0]:
            raise ValueError(f"{name}: independent End(C^2) character {direct}, branching gives {totals[0]}")
        output[name] = dict(weight=weight, length=case["length"], parity=case["colour_parity"],
             branches=len(case["active"]), cochains=direct, classical_dimension_upper=totals[1],
             cartan_kernel_upper=totals[2], remaining=remaining)
    return dict(prime=prime, sectors=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path)
    parser.add_argument("--raw", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--prime", type=int, choices=(1000003, 1000033, 2147483647), required=True)
    parser.add_argument("--case", action="append", choices=tuple(case[0] for case in CASES))
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--additional-records", type=Path, action="append", default=[],
                        help="merge disjoint or identical active-sector records before replay")
    parser.add_argument("--resolve", type=Path)
    args = parser.parse_args()
    if (not __debug__ or bool(args.raw) != bool(args.engine) or not (args.engine or args.resolve)
            or (args.additional_records and not args.resolve)):
        parser.error("use --engine with --raw or replay with --resolve")
    cases = {name: row for name, row in plan().items() if not args.case or name in args.case}
    if args.engine:
        records = compute(cases, args.engine, args.raw, args.workers, args.prime)
        gzip_json(args.records, records)
    if args.resolve:
        records = json.loads(gzip.decompress(args.records.read_bytes()))
        for path in args.additional_records:
            extra = json.loads(gzip.decompress(path.read_bytes()))
            if any(records[k] != extra[k] for k in records.keys() & extra.keys()):
                raise ValueError(f"conflicting active-sector records in {path}")
            records.update(extra)
        result = resolve(cases, records, args.prime)
        args.resolve.write_text(json.dumps(result, indent=2) + "\n")
