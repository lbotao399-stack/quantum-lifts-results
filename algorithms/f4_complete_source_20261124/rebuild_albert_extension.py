"""Rebuild only the literal Albert boundary/product columns, with no colour selection map."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
from hashlib import sha256
import json
from multiprocessing import get_context
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

from albert_tensor_runner import AlbertRunner


def initialize(inputs, library, prime):
    global runner
    runner = AlbertRunner(Path(inputs), Path(library), prime)


def evaluate(item):
    position, point = item
    started = process_time()
    return position, runner.evaluate(point), process_time() - started


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--prime", type=int, choices=(1009, 1013), default=1009)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--reuse-prefix", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.workers <= 2:
        parser.error("one or two workers are allowed")
    global_directory = args.inputs.parent / "global_boundaries"
    old = json.loads((global_directory / "channels-3060.json").read_text())
    mixed = json.loads((global_directory / "mixed-minor.json").read_text())
    fixed = json.loads((args.inputs / f"sixth-product-{args.prime}.json").read_text())
    points = old["rows"] + mixed["new_rows"] + fixed["new_points"]
    assert len(points) == len(set(points)) == 3017
    value_path = Path(str(args.out) + ".npy")
    done_path = Path(str(args.out) + ".complete.npy")
    receipt_path = Path(str(args.out) + ".json")
    hashes = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
              for name in ("albert_tensor_runner.cpp", "albert_tensor_runner.py")}
    if receipt_path.exists():
        state = json.loads(receipt_path.read_text())
        assert state["points"] == points and state["implementation_sha256"] == hashes
        values = np.lib.format.open_memmap(value_path, mode="r+")
        done = np.lib.format.open_memmap(done_path, mode="r+")
    else:
        values = np.lib.format.open_memmap(value_path, mode="w+", dtype="<u2", shape=(3017, 50))
        done = np.lib.format.open_memmap(done_path, mode="w+", dtype=np.uint8, shape=(3017,))
        done[:] = 0
        state = {"schema": "f4-literal-albert-extension-v1", "prime": args.prime,
                 "points": points, "shape": [3017, 50], "complete": False,
                 "implementation_sha256": hashes, "worker_cpu_seconds_this_stage": 0,
                 "reused_literal_tensor_rows": 0}
        if args.reuse_prefix:
            prior = json.loads(Path(str(args.reuse_prefix) + ".json").read_text())
            assert prior["prime"] == args.prime and prior["points"] == points
            assert prior["hashes"]["tensor_runner_source"] == hashes["albert_tensor_runner.cpp"]
            assert prior["hashes"]["tensor_runner_python"] == hashes["albert_tensor_runner.py"]
            prior_done = np.load(Path(str(args.reuse_prefix) + ".complete.npy"))
            prior_values = np.load(Path(str(args.reuse_prefix) + ".extension.npy"), mmap_mode="r")
            values[prior_done.astype(bool)] = prior_values[prior_done.astype(bool)]
            done[:] = prior_done
            state["reused_literal_tensor_rows"] = int(done.sum())
            state["reuse_scope"] = "Only the literal AlbertRunner output columns were copied. No baseline matrix values or colour-selection operations enter this standalone rebuild."
        values.flush(); done.flush()
        receipt_path.write_text(json.dumps(state, indent=2) + "\n")
    tasks = [(i, point) for i, point in enumerate(points) if not done[i]]
    started = last = perf_counter()
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork"), initializer=initialize,
                             initargs=(str(args.inputs), str(args.library), args.prime)) as pool:
        for position, value, cpu in pool.map(evaluate, tasks, chunksize=2):
            values[position] = value; done[position] = 1
            state["worker_cpu_seconds_this_stage"] += cpu
            if perf_counter() - last > 30:
                values.flush(); done.flush()
                state["completed_rows"] = int(done.sum())
                receipt_path.write_text(json.dumps(state, indent=2) + "\n")
                print(json.dumps({"completed_rows": int(done.sum()), "total_rows": 3017}), flush=True)
                last = perf_counter()
    values.flush(); done.flush()
    assert done.all()
    raw = values.tobytes()
    compressed = Path(str(args.out) + ".u16.gz")
    compressed.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))
    state.update(complete=True, completed_rows=3017, native_entries=3017 * 50,
                 residue_u16le_sha256=sha256(raw).hexdigest(),
                 gzip_sha256=sha256(compressed.read_bytes()).hexdigest(),
                 scope="All 49 literal Albert boundaries and the closed product are freshly evaluated at the 3017 published points. This does not independently reconstruct the old 3002 boundary block or certify the full 3017 minor.")
    receipt_path.write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps({k: v for k, v in state.items() if k != "points"}), flush=True)


if __name__ == "__main__":
    main()
