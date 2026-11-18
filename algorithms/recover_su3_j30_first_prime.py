"""Recover the SU(3) level-30 joint solution from 50 saved first-prime stages.

The zero-free-coordinate solution is checkpointed before multiplying the original
map blocks. Each residual is recomputed from the block hashes in the stage that
consumed it, so elimination and residual multiplication can be resumed separately.
"""
import argparse
import hashlib
import io
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np


PRIME = 1000003
FIRST = 28311
JOINT = 77639


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            result.update(block)
    return result.hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    if isinstance(value, np.ndarray):
        with temporary.open("wb") as stream:
            np.save(stream, value, allow_pickle=False)
    else:
        temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def stages(root):
    previous = []
    used = set()
    records = []
    first = root / "su3_coupled_completion_v1_20261117/solve/calculation"
    second = root / "su3_second_solve_v1_20261118/calculation"
    for order, parent, size, variables in ((1, first, 13, FIRST), (2, second, 37, JOINT)):
        folders = sorted(parent.glob("stage*_target*"))
        assert len(folders) == size
        for index, folder in enumerate(folders):
            path = folder / "result.json"
            record = json.loads(path.read_text())
            manifest = record["manifest"]
            assert folder.name == "stage%02d_target%02d" % (index, manifest["target"])
            assert manifest["previous_stages"] == previous
            assert manifest["prime"] == PRIME and manifest["variables"] == variables
            assert record["consistent"] and len(record["pivots"]) == record["rank"]
            pivots = set(record["pivots"])
            assert len(pivots) == record["rank"] and used.isdisjoint(pivots)
            assert all(0 <= p < variables for p in pivots)
            used.update(pivots)
            previous.append(digest(path))
            records.append((order, folder, record))
        assert len(used) == (18003 if order == 1 else 50147)
    return previous, records, sorted(set(range(JOINT)) - used)


def recover(records, free):
    solution = np.zeros(JOINT, dtype=np.int64)
    for order, folder, record in reversed(records):
        for name in ("rows", "columns"):
            assert digest(folder / (name + ".npy")) == record[name + "_sha256"]
        columns = np.load(folder / "columns.npy", allow_pickle=False)
        rows = np.load(folder / "rows.npy", allow_pickle=False, mmap_mode="r")
        assert columns[-1] == (FIRST if order == 1 else JOINT)
        assert np.all(np.diff(columns) > 0) and columns[0] >= 0
        assert rows.shape == (record["rank"], len(columns))
        assert len(columns) * (PRIME - 1) ** 2 < 2 ** 63
        known = solution[columns[:-1]]
        for start in range(0, len(rows), 128):
            block = rows[start:start+128].astype(np.int64)
            assert np.all((block >= 0) & (block < PRIME))
            pivots = record["pivots"][start:start + len(block)]
            solution[pivots] = (block[:, -1] - block[:, :-1] @ known) % PRIME
        known = solution[columns[:-1]]
        for start in range(0, len(rows), 128):
            block = rows[start:start+128].astype(np.int64)
            assert not np.any((block[:, :-1] @ known - block[:, -1]) % PRIME)
    assert len(free) == 27492 and not np.any(solution[free])
    return solution


def checked_array(path, expected):
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == expected, path
    values = np.load(io.BytesIO(data), allow_pickle=False)
    assert values.dtype.kind in "iu" and np.all((values >= 0) & (values < PRIME))
    return values


def replay(root, runtime, records, solution, output, workers, previous=None):
    first_maps = root / "su3_coupled_map_v2_20261117"
    roots = {(1, "Q0"): first_maps, (2, "Q0"): root / "su3_second_maps_v1_20261118",
             (2, "Q1"): root / "su3_second_Q1_v1_20261118"}
    plans = {1: "su3_correction_profiles_v1.json", 2: "su3_second_correction_profiles_v1.json"}
    offsets = {}
    dimensions = {}
    for order, name in plans.items():
        cursor = 0
        for source in json.loads((runtime / name).read_text())["sources"]:
            offsets[order, source["index"]] = cursor
            dimensions[order, source["index"]] = source["dimension"]
            cursor += source["dimension"]
        assert cursor == (FIRST if order == 1 else JOINT - FIRST)
    parts = {}
    for (order, kind), folder in roots.items():
        for marker in (folder / "maps").glob("profile*/**/manifest.json"):
            manifest = json.loads(marker.read_text())
            assert manifest["prime"] == PRIME
            key = order, kind, manifest["source"]["index"]
            parts.setdefault(key, []).append((manifest["start_row"], manifest["stop_row"], marker.parent))

    def multiply(item):
        order, stage, record = item
        manifest = record["manifest"]
        target = manifest["target"]
        width = manifest["columns"]
        name = "residual_%d_%02d.json" % (order, target)
        if previous and (previous / name).is_file():
            receipt = json.loads((previous / name).read_text())
            assert receipt == dict(order=order, target=target, coordinates=width,
                                   blocks=len(manifest["blocks"]), nonzero=0,
                                   stage_sha256=digest(stage / "result.json"))
            save(output / name, receipt)
            return receipt
        total = np.zeros(width, dtype=np.int64)
        rhs_root = first_maps if order == 1 else root / "su3_second_rhs_v1_20261118"
        rhs = checked_array(rhs_root / "rhs" / ("profile%02d" % target) / "rhs.npy", manifest["rhs_sha256"])
        assert rhs.shape == total.shape
        coverage = {}
        for block in manifest["blocks"]:
            source, first, last = block["source"], block["start"], block["stop"]
            kind = block.get("kind", "Q0")
            source_order = 2 if order == 2 and kind == "Q0" else 1
            offset = offsets[source_order, source] + (FIRST if source_order == 2 else 0)
            if order == 2:
                assert block["global_start"] == offset + first
            matches = [p for a, b, p in parts[order, kind, source] if a <= first < last <= b]
            assert len(matches) == 1
            path = matches[0] / ("t%02d_%05d_%05d.npy" % (target, first, last))
            values = checked_array(path, block["sha256"])
            assert values.shape == (last - first, width)
            assert (last - first) * (PRIME - 1) ** 2 + PRIME < 2 ** 63
            total = (total + solution[offset + first:offset + last] @ values) % PRIME
            coverage.setdefault((source_order, source), []).append((first, last))
        for key, spans in coverage.items():
            cursor = 0
            for first, last in sorted(spans):
                assert first == cursor
                cursor = last
            assert cursor == dimensions[key]
        assert not np.any((total - rhs) % PRIME), (order, target)
        receipt = dict(order=order, target=target, coordinates=width, blocks=len(manifest["blocks"]),
                       nonzero=0, stage_sha256=digest(stage / "result.json"))
        save(output / name, receipt)
        return receipt

    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(multiply, records))
    assert [sum(r["coordinates"] for r in results if r["order"] == n) for n in (1, 2)] == [34901, 65406]
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--replay-from", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    hashes, records, free = stages(args.root)
    if args.replay_from:
        checkpoint = json.loads((args.replay_from / "recovery.json").read_text())
        assert checkpoint["prime"] == PRIME and checkpoint["stages"] == hashes
        assert checkpoint["variables"] == JOINT and checkpoint["rank"] == 50147
        solution_path = args.replay_from / "particular.npy"
        free_path = args.replay_from / "free_variables.npy"
        assert digest(solution_path) == checkpoint["solution_sha256"]
        assert digest(free_path) == checkpoint["free_sha256"]
        assert np.array_equal(np.load(free_path, allow_pickle=False), free)
        solution = np.load(solution_path, allow_pickle=False).astype(np.int64)
        assert solution.shape == (JOINT,) and np.all((solution >= 0) & (solution < PRIME))
        assert not np.any(solution[free])
    else:
        solution = recover(records, free)
    save(args.output / "particular.npy", solution.astype(np.uint32))
    save(args.output / "free_variables.npy", np.asarray(free, dtype=np.int64))
    record = dict(prime=PRIME, stages=hashes, variables=JOINT, rank=50147, free=len(free),
                  job=os.environ.get("SLURM_JOB_ID"), solution_sha256=digest(args.output / "particular.npy"),
                  free_sha256=digest(args.output / "free_variables.npy"))
    save(args.output / "recovery.json", record)
    runtime = args.root / "package/work/sp3-native.nosync"
    residuals = replay(args.root, runtime, records, solution, args.output, args.workers, args.replay_from)
    save(args.output / "result.json", dict(record, residuals=residuals))


if __name__ == "__main__":
    main()
