"""Classical Cartan-kernel minors for connected Levi factors of the F4 source."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import product
import json
import os
from pathlib import Path
import subprocess


SHAPE = (4, 4, 4, 2)
PRIME_SEEDS = {1000003: 372071, 1000033: 371993}


def branches(central_rank):
    """Remove at most one horizontal strip for each central scalar colour."""
    shapes = {SHAPE}
    for _ in range(central_rank):
        shorter = set()
        for shape in shapes:
            choices = [range(low, high + 1) for high, low in zip(shape, (*shape[1:], 0))]
            shorter.update(tuple(k for k in row if k) for row in product(*choices))
        shapes = shorter
    return sorted(shapes, key=lambda shape: (sum(shape), shape))


def taylor_weight(shape):
    # Three ordinary first-scalar directions followed by the two odd first gauginos.
    scalars = (shape + (0, 0, 0))[:3]
    dotted = tuple(sum(row >= j for row in shape[3:]) for j in (1, 2))
    assert sum(scalars + dotted) == sum(shape)
    return scalars + dotted


def validate(entry, group, shape, prime, weight=None):
    weight = taylor_weight(shape) if weight is None else tuple(weight)
    assert entry["gauge"] == group and tuple(entry["sector"]) == weight
    assert entry["length"] == sum(shape) - 1 and entry["seed"] == PRIME_SEEDS[prime]
    assert entry["prime"] == prime and entry["mode"] == "classical_Cartan_kernel"
    dimension = entry["cochains"]
    if not dimension:
        return 0, 0, 0, 0, 0
    for key, height, width in (
        ("incoming", dimension, entry["parent"]),
        ("outgoing", entry["next"], dimension),
        ("outgoing_and_Cartan", entry["stacked_rows"], dimension),
    ):
        minor = entry[key]
        rank = minor["rank"]
        assert 0 < minor["det"] < prime
        assert len(set(minor["rows"])) == len(minor["rows"]) == rank
        assert len(set(minor["columns"])) == len(minor["columns"]) == rank
        assert all(0 <= row < height for row in minor["rows"])
        assert all(0 <= column < width for column in minor["columns"])
    incoming = entry["incoming"]["rank"]
    outgoing = entry["outgoing"]["rank"]
    stacked = entry["outgoing_and_Cartan"]["rank"]
    bound = dimension - incoming - stacked
    assert bound >= 0 and stacked >= outgoing
    return dimension, incoming, outgoing, stacked, bound


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", type=Path, help="executable built from su2_cartan_kernel.cpp")
    parser.add_argument("directory", type=Path)
    parser.add_argument("--groups", nargs="+", choices=("SU3", "Sp2", "SO7", "Sp3"),
                        default=("SU3", "Sp2"))
    parser.add_argument("--primes", nargs="+", type=int, choices=tuple(PRIME_SEEDS), default=tuple(PRIME_SEEDS))
    parser.add_argument("--min-rows", type=int, default=1)
    parser.add_argument("--max-rows", type=int, default=4)
    parser.add_argument("--min-degree", type=int, default=1)
    parser.add_argument("--max-degree", type=int, default=14)
    parser.add_argument("--shape", type=int, nargs="+", help="one ordinary Levi diagram")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    if not __debug__ or min(args.workers, args.threads) < 1:
        parser.error("enable assertions and use positive parallelism")
    args.directory.mkdir(parents=True, exist_ok=True)
    executable = args.engine.resolve()
    environment = dict(os.environ, PAR_THREADS=str(args.threads))
    cases = [(group, shape, prime) for group in args.groups
             for shape in branches(2 if group in ("SU3", "Sp2") else 1)
             if args.min_rows <= len(shape) <= args.max_rows
             and args.min_degree <= sum(shape) <= args.max_degree
             and (args.shape is None or shape == tuple(args.shape))
             for prime in args.primes]
    if not cases:
        parser.error("no diagram matches the requested group and filters")

    def check(case):
        group, shape, prime = case
        name = "-".join(map(str, (group, *shape, prime)))
        target = args.directory / (name + ".json")
        if target.exists():
            entry = json.loads(target.read_bytes())
        else:
            command = [str(executable), *map(str, taylor_weight(shape)), str(sum(shape) - 1),
                       str(prime), str(PRIME_SEEDS[prime]), group]
            entry = json.loads(subprocess.check_output(command, env=environment))
            validate(entry, group, shape, prime)
            temporary = target.with_suffix(f".{os.getpid()}.part")
            temporary.write_text(json.dumps(entry, separators=(",", ":")) + "\n")
            temporary.replace(target)
        ranks = validate(entry, group, shape, prime)
        print(name, *ranks, flush=True)

    with ThreadPoolExecutor(args.workers) as pool:
        for _ in pool.map(check, cases):
            pass


if __name__ == "__main__":
    main()
