"""Contract the two-exchange classical columns and the terminal graphs for completed Ward reductions."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import hashlib
from itertools import combinations
import json
from multiprocessing import get_context
import os
from pathlib import Path
import sys
import tempfile
from zipfile import ZipFile

import numpy as np
import integer_section as integer
from other_crossing import ROOT, FILES, EXPECTED, COLOURS, insertion_sign

PAIRS = tuple(combinations(range(4), 2))


def point(i):
    with np.errstate(all="ignore"):
        x = (integer.sample(i) @ integer.native26().reshape(52, 676)).reshape(4, 26, 26)
        low, high, exchange, scales = integer.tensors()
        assert scales == (28, 252, 84)
        left = -integer.halves(
            x.transpose(0, 2, 1), low, [integer.CROSSED[:2]]
        ).reshape((26,) * 4)
        right = integer.halves(x, high, [integer.CROSSED[2:]]).reshape((26,) * 4)
        m, mt = integer.action(exchange), integer.action(exchange.T)
        left = np.stack([integer.exact_apply(mt, left, p).ravel() for p in PAIRS])
        right = np.stack([integer.exact_apply(m, right, p).ravel() for p in PAIRS])
        # The first recorded exchange acts on the right column first.
        products = integer.gram(left, right).T
    return [i, list(map(int, products.ravel()))]


def graph_key(g):
    return json.dumps(
        [g[k] for k in ("E", "d", "cp", "cm", "cubics")], separators=(",", ":")
    )


def poles(task):
    import native

    idx, g = task
    key = tuple(map(tuple, g["E"])), tuple(g["d"]), tuple(g["cp"]), tuple(g["cm"])
    _, owners, tree = native.plan(*key)
    sign = insertion_sign(g, owners, tree)
    return [idx, [str(sign * native.run(g, j)) for j in (0, 1)]]


def selected(banks):
    keys = {}
    for item in banks:
        if any(v["ordinary"] for v in item["banks"].values()):
            continue
        for data in item["banks"].values():
            for g in data["banks"]:
                keys.setdefault(graph_key(g), g)
    return dict(sorted(keys.items()))


def run(args):
    if args.part == "columns":
        with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
            out = list(pool.map(point, range(720), chunksize=2))
        args.output.write_bytes(
            gzip.compress(json.dumps(out, separators=(",", ":")).encode(), mtime=0)
        )
        return
    banks = json.loads(gzip.decompress(args.banks.read_bytes()))
    if isinstance(banks, dict):
        banks = banks["kernel"]
    keys = selected(banks)
    todo = [
        (i, g) for i, g in enumerate(keys.values()) if i % args.shards == args.shard
    ]
    os.environ["P"] = "0"
    with tempfile.TemporaryDirectory() as temp, ZipFile(args.packet) as z:
        root = Path(temp) / ROOT
        for name, digest in zip(FILES, EXPECTED, strict=True):
            raw = z.read(ROOT + "/" + name)
            assert hashlib.sha256(raw).hexdigest() == digest
            f = root / name
            f.parent.mkdir(exist_ok=True, parents=True)
            f.write_bytes(raw)
        sys.path.insert(0, str(root / COLOURS))
        import native

        assert native.P == 0
        saved = (
            dict(json.loads(args.output.read_bytes())) if args.output.exists() else {}
        )
        with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
            todo = [pool.submit(poles, t) for t in todo if t[0] not in saved]
            for f in as_completed(todo):
                idx, value = f.result()
                saved[idx] = value
                tmp = args.output.with_suffix(".part")
                tmp.write_text(
                    json.dumps(sorted(saved.items()), separators=(",", ":")) + "\n"
                )
                tmp.replace(args.output)
                if len(saved) % 10 == 0:
                    print("native rational graph count", len(saved), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("part", choices=["columns", "poles"])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--banks", type=Path)
    ap.add_argument("--packet", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    args = ap.parse_args()
    assert 0 <= args.shard < args.shards
    run(args)
