"""Compute the charge-conjugation ranks at derivative grade two for incoming F1 partners."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


PRIME = 2147483647


def save(path, record):
    pending = path.with_suffix(".part")
    pending.write_text(json.dumps(record, indent=2) + "\n")
    pending.replace(path)


def calculate(engine, directory, threads, job):
    directory.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for level, length, tail in [(30, 7, 3), (32, 8, 4)]:
        for index in reversed(range(5)):
            weight = [tail + 2 + index, 0, 0, 4, 4 - index]
            seed = 611209 + level * 61 + index * 1009
            path = directory / f"j{level}_sector_{index}.json"
            if path.exists():
                row = json.loads(path.read_bytes())
            else:
                args = [str(engine.resolve()), "SU3", *map(str, weight), str(length), str(seed), str(PRIME)]
                row = json.loads(subprocess.check_output(args, env=dict(os.environ, PAR_THREADS=str(threads))))
            if (row["mode"], row["gauge"], row["sector"], row["length"], row["seed"], row["prime"]) != (
                    "charge_conjugation", "SU3", weight, length, seed, PRIME):
                raise ValueError(f"Unexpected native ranks for {path.name}")
            ranks = [row[k] for k in ("even_in", "odd_in", "even_out", "odd_out")]
            if min(ranks) < 0 or sum(ranks) > row["cochains"]:
                raise ValueError(f"Incompatible adjacent maps in {path.name}")
            save(path, row)
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            save(directory / "partial.json", dict(job=job, files=hashes))
            print(level, index, row["cochains"], ranks, flush=True)
    (directory / "partial.json").replace(directory / "receipt.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--job", default="local")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("threads must be positive")
    calculate(args.engine, args.output, args.threads, args.job)
