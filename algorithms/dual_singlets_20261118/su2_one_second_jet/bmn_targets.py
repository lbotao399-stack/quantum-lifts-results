"""Test the spin-singlet targets of the SU(2) BMN second obstruction through level 30."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess

PRIMES_AND_SEEDS = ((1000003, 772801), (1000033, 843119))


def weights():
    for total in range(10):
        for first in range(total + 1):
            for second in range(min(total - first, first) + 1):
                last = total - first - second
                if 0 <= last <= second:
                    yield (first, second, last)


def record(weight, results):
    total = sum(weight)
    columns = ("Y", "dim", "rank_in", "rank_out", "rho")
    values = []
    for result, (prime, seed) in zip(results, PRIMES_AND_SEEDS, strict=True):
        assert result["gauge"] == "SU2" and result["sector"] == [*weight, 2, 2]
        assert (result["prime"], result["seed"], result["j"]) == (prime, seed, 2 * total + 12)
        block = result["blocks"]
        assert [v["Y"] for v in block] == [total + 2, total + 3]
        for v in block:
            assert v["dim"] == v["rank_in"] + v["rank_out"] + v.get("rho", 0)
            assert "rho" in v or v["dim"] == 0
        assert block[1].get("rho", 0) == 0
        values.append([[v.get(k, 0) for k in columns] for v in block])
    assert values[0] == values[1]
    return dict(flavour=list(weight), integer_dimension_and_modular_ranks=values[0])


def run(args):
    def one(weight):
        total = sum(weight)
        results = []
        for prime, seed in PRIMES_AND_SEEDS:
            cmd = [str(args.engine.resolve()), "chain", "SU2", *map(str, weight),
                   "2", "2", str(total + 2), str(total + 3), str(seed), str(prime)]
            env = dict(os.environ, PAR_THREADS="1", A_CHECK="1", A_VERB="0")
            env.pop("A_NORHO", None)
            results.append(json.loads(subprocess.check_output(cmd, env=env)))
        return record(weight, results)

    with ThreadPoolExecutor(args.workers) as pool:
        rows = list(pool.map(one, weights()))
    return dict(primes_and_seeds=PRIMES_AND_SEEDS,
                column_order=["letter_number", "cochains", "incoming_rank", "outgoing_rank", "cartan_rank_on_cycles"],
                dot_degrees=[2, 2], levels=[12, 30], records=rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("engine", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if not __debug__ or args.workers < 1:
        parser.error("Enable assertions and supply a positive worker count")
    print(json.dumps(run(args), indent=2))
