"""Compute adjacent classical ranks at one lower USp(4) BMN source."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

try:
    from b2_schur_complex import multiplicity, transpose
except ModuleNotFoundError as error:
    if error.name != "b2_schur_complex":
        raise
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from b2_schur_complex import multiplicity, transpose


def calculate(engine, weight, length, prime, seed, threads):
    shape = transpose(tuple(weight))
    dimension = multiplicity(shape, length)
    args = [str(engine.resolve()), "chain", "Sp2", *map(str, weight), "0", "0",
            str(length), str(length), str(seed), str(prime)]
    env = dict(os.environ, PAR_THREADS=str(threads), A_NORHO="1", A_VERB="1")
    record = json.loads(subprocess.check_output(args, env=env))
    block, = record["blocks"]
    if (record["gauge"], record["sector"], record["prime"], record["seed"], block["Y"], block["dim"]) != (
            "Sp2", [*weight,0,0], prime, seed, length, dimension):
        raise ValueError("native and independent character calculations disagree")
    if not dimension-block["rank_in"]-block["rank_out"] == block["h"] >= 0:
        raise ValueError("adjacent ranks violate the cochain dimension")
    return dict(group="USp(4)", weight=weight, level=record["j"], length=length, prime=prime,
                seed=seed, cochains=dimension, incoming=block["rank_in"], outgoing=block["rank_out"],
                rational_cohomology_upper=block["h"])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--weight", type=int, nargs=3, required=True)
    p.add_argument("--length", type=int, required=True)
    p.add_argument("--seed", type=int, default=938721)
    p.add_argument("--prime", type=int, choices=(1000003,1000033,2147483647), default=2147483647)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    r = calculate(a.engine, a.weight, a.length, a.prime, a.seed, a.threads)
    a.out.write_text(json.dumps(r, indent=2)+"\n")
