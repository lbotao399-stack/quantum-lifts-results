"""Evaluate both marks and both ordered exchanges in the pinned eight-Albert packet."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from hashlib import sha256
from itertools import chain
import json
from multiprocessing import get_context
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
import integer_section as integer

HERE = Path(__file__).resolve().parent
ROOT = "f4-j36-eight-Albert-certificate"
COLOURS = "f4-j36-eight-cubic"
POLES = (-4 * 2392 * 7113106, 2 * 2299 * 1850212)
NORMALIZATION = 2 * Fraction(28, 3) ** 4
FILES = (COLOURS + "/native.py", COLOURS + "/bank-A-reduced.json",
         "fixtures/color-basis-native.json", "f4-j36-nested1-seal/cubic.json")
PACKET_SHA = "b668a97232466286970c3c3b7df49128d5b98e8e0a4b00e23f57fa75563f7972"


def column(point, dense=False):
    x = (integer.sample(point, dense) @ integer.native26().reshape(52, 676)).reshape(4, 26, 26)
    low, high, exchange, scales = integer.tensors()
    left_arm, right_arm = integer.CROSSED[:2], integer.CROSSED[2:]
    left = -integer.halves(x.transpose(0, 2, 1), low, [left_arm])
    right = integer.halves(x, high, [right_arm]).reshape((26,) * 4)
    transfer = integer.action(exchange)
    results = []
    for axes in [((0, 1), (1, 2)), ((1, 2), (0, 1))]:
        value = right
        for pair in axes:
            value = integer.exact_apply(transfer, value, pair)
        results.append(int(integer.gram(left, value.reshape(1, -1))[0, 0]))
    assert scales == (28, 252, 84)
    return [point, *results]


def insertion_sign(graph, owners, tree):
    edges = sorted(set(graph["cp"]) | set(graph["cm"]))
    owners = {e: graph["E"][e][end] for e, end in zip(edges, owners, strict=True)}

    def visit(node):
        if isinstance(node, int):
            return [e for e in edges if owners[e] == node]
        return list(chain.from_iterable(visit(child) for child in node))

    order = visit(tree)
    assert sorted(order) == edges
    inversions = 0
    for occupied in (graph["cp"], graph["cm"]):
        positions = [occupied.index(e) for e in order if e in occupied]
        inversions += sum(x > y for i, x in enumerate(positions) for y in positions[i + 1:])
    return (-1) ** inversions


def pole_row(task):
    import native
    refs, graph = task
    args = tuple(map(tuple, graph["E"])), tuple(graph["d"]), tuple(graph["cp"]), tuple(graph["cm"])
    _, owners, tree = native.plan(*args)
    poles = [insertion_sign(graph, owners, tree) * native.run(graph, j) for j in (0, 1)]
    return dict(entries=refs, ordered_poles=[str(v) for v in poles], weight=graph["w"])


def run_poles(folder, workers, output):
    os.environ["P"] = "0"
    sys.path.insert(0, str(folder / COLOURS))
    import native
    assert native.P == 0
    banks = json.loads((folder / COLOURS / "bank-A-reduced.json").read_bytes())
    unique = {}
    for mark in ("left", "right"):
        assert not banks[mark]["ordinary"]
        for i, graph in enumerate(banks[mark]["banks"]):
            key = json.dumps(graph, sort_keys=True)
            unique.setdefault(key, ([], graph))[0].append([mark, i])
    with ProcessPoolExecutor(workers, mp_context=get_context("fork")) as pool:
        futures = [pool.submit(pole_row, t) for t in unique.values()]
        from concurrent.futures import as_completed
        rows = []
        for f in as_completed(futures):
            rows.append(f.result())
            print("rational graph", rows[-1]["entries"], rows[-1]["ordered_poles"], flush=True)
            output.write_text(json.dumps(rows, indent=2) + "\n")


def run(args):
    if args.mode == "columns":
        integer.tensors()
        with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
            rows = list(pool.map(column, range(720), chunksize=2))
        args.output.write_text(json.dumps(rows, separators=(",", ":")) + "\n")
        return
    assert args.packet
    with TemporaryDirectory() as temp, ZipFile(args.packet) as source:
        root = Path(temp) / ROOT
        if sha256(args.packet.read_bytes()).hexdigest() != PACKET_SHA:
            # A reduced transport input must contain precisely the four pinned members.
            assert {s for s in source.namelist()} == {ROOT + "/" + s for s in FILES}
        for name, digest in zip(FILES, EXPECTED, strict=True):
            raw = source.read(ROOT + "/" + name)
            assert sha256(raw).hexdigest() == digest
            out = root / name
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(raw)
        run_poles(root, args.workers, args.output)


EXPECTED = (
    "96d6cde098b6281d0dad13951d9603b4f9d0f86f9117226b1efab9f92dd0648b",
    "1b1d298cf63b04e4b95239758af57564c3e42de6c1682dcc8cb2e7d0dc1da502",
    "030c2bf73e82a33b59721dea76bd3478d24ee7e01fc6bf805a6c97d0ba29e71c",
    "03ad3705cc69c129c53b4c6b2e43891e910b2c619a7f6a9bf24070b821937221",
)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=("columns", "poles"))
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--workers", default=1, type=int)
    ap.add_argument("--packet", type=Path)
    run(ap.parse_args())
