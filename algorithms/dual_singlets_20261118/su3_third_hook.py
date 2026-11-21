"""Compute the level-36 SU(3) third-order Taylor target one classical map at a time."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from taylor_hook_cohomology import saved, weights
from trace_descent import dimensions
from two_flavour_hooks import schur_colour_dimensions


PARTITION = (6, 6, 1, 1, 1)
SECTORS = weights(6, 3, False)


def cochains(sector, length):
    degree = (*sector[3:], *sector[:3])
    return [sum(dimensions(degree, length, 3, colour).values()) for colour in ("SU+", "SU-")]


def check_image(record, sector, length, prime):
    if (record["mode"], record["gauge"], record["sector"], record["length"], record["prime"]) != (
            "complete_colour_image", "SU3", list(sector), length, prime):
        raise ValueError("image has a different sector or frame type")
    source = cochains(sector, length)
    target = cochains(sector, length + 1)
    if (sum(source), sum(target)) != (record["cochains"], record["target_evaluations"]):
        raise ValueError("incomplete source or target in the finite-colour character")
    ranks = []
    for i, parity in enumerate(("even", "odd")):
        minor = record[parity]
        rank = minor["rank"]
        if not 0 <= rank <= min(source[i], target[i]) or not 0 < minor["det"] < prime:
            raise ValueError("invalid rank or determinant")
        for key, count in [("columns", sum(source)), ("rows", sum(target))]:
            values = minor[key]
            if len(values) != rank or len(set(values)) != rank or any(k not in range(count) for k in values):
                raise ValueError("invalid minor coordinates")
        ranks.append(rank)
    return ranks


def map_path(folder, length, index):
    return folder / f"image_{length}_{index}.json"


def compute(executable, folder, lengths, threads, prime, seed):
    env = dict(os.environ, PAR_THREADS=str(threads))
    folder.mkdir(exist_ok=True, parents=True)
    for length in lengths:
        for i, sector in enumerate(SECTORS):
            path = map_path(folder, length, i)
            map_seed = seed + 1009 * i + 7919 * length
            if path.is_file():
                image = json.loads(path.read_bytes())
            else:
                args = [str(executable.resolve()), "SU3", *map(str, sector), str(length),
                        str(map_seed), str(prime), "complete-image"]
                image = json.loads(subprocess.check_output(args, env=env))
            if image["seed"] != map_seed:
                raise ValueError("seed mismatch in saved image")
            ranks = check_image(image, sector, length, prime)
            saved(path, image)
            print(length, i, *ranks, flush=True)


def resolve(folder, lengths, prime):
    blocks = []
    for length in sorted(set(lengths)):
        cochain = list(schur_colour_dimensions(PARTITION, length, 3))
        character = [0, 0]
        incoming, outgoing = [0, 0], [0, 0]
        for i, sector in enumerate(SECTORS):
            sign = (-1)**i
            dim = cochains(sector, length)
            for k in range(2):
                character[k] += sign * dim[k]
            for y, result in [(length-1, incoming), (length, outgoing)]:
                if y == sum(PARTITION):
                    continue
                raw = json.loads(map_path(folder, y, i).read_bytes())
                ranks = check_image(raw, sector, y, prime)
                for k in range(2):
                    result[k] += sign * ranks[k]
        h = [dim-ri-ro for dim, ri, ro in zip(cochain, incoming, outgoing)]
        if character != cochain or min(incoming + outgoing + h) < 0:
            raise ValueError("Schur subtraction and independent character disagree")
        blocks.append(dict(target_length=length, source_length=length+2, cochains=cochain,
                           incoming=incoming, outgoing=outgoing, cohomology_mod_prime=h))
    return dict(group="SU(3)", source_weight=[8, 5, 5], order=3, partition=PARTITION,
                prime=prime, parities=["even", "odd"], blocks=blocks)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder", type=Path, required=True)
    p.add_argument("--engine", type=Path)
    p.add_argument("--maps", nargs="+", type=int, choices=range(1, 15))
    p.add_argument("--resolve", nargs="+", type=int, choices=range(2, 16))
    p.add_argument("--out", type=Path)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--prime", type=int, choices=(1000003, 1000033, 2147483647), default=2147483647)
    p.add_argument("--seed", type=int, default=913873)
    a = p.parse_args()
    if not (a.maps or a.resolve) or bool(a.maps) != bool(a.engine) or bool(a.resolve) != bool(a.out):
        p.error("supply --maps with --engine or --resolve with --out")
    if a.maps:
        compute(a.engine, a.folder, a.maps, a.threads, a.prime, a.seed)
    if a.resolve:
        saved(a.out, resolve(a.folder, a.resolve, a.prime))
