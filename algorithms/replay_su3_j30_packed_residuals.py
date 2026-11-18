"""Independently multiply the level-30 saved solution by all three packed SU(3) maps."""

import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np


PRIME = 1000033
PLANS = {
    "first": ("su3_correction_profiles_v1.json", "b4eb0388bbaa1a6dc655e4b2d877ba48e771d00a89cc7a76d48f4278219570eb"),
    "second_Q1": ("su3_second_quantum_profiles_v1.json", "e0a9ee4181b3c12ae0f0b0bdb0367fd30c4d8d55c535ba588345692966adb462"),
    "second_Q0": ("su3_second_correction_profiles_v1.json", "b8c669fa244e851f38a25b91186c3091e98a939935242805b0f71a9f8c00d609"),
}


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_text())


def vector(path, digest, size):
    assert sha(path) == digest
    v = np.load(path, allow_pickle=False)
    assert v.shape == (size,) and v.dtype == np.dtype("<u4") and np.all(v < PRIME)
    return v.astype(np.int64)


def replay(root, workers):
    runtime = root / "runtime/work/sp3-native.nosync"
    output = root / "output"
    joint = vector(output / "second_solve_compact/particular.npy",
                   "a3f3ea9e9448299610f5313c7a4271184b099850b7696989c5bf8f5ea5dd6f39", 77639)
    first = vector(output / "first_solve_compact/particular.npy",
                   "7471532cfc2ba11e79f524435625876376d8b91290326d87567ce43bd2d2ff33", 28311)
    plans, offsets = {}, {}
    for kind, (name, digest) in PLANS.items():
        path = runtime / name
        assert sha(path) == digest
        plan = plans[kind] = read(path)
        offset = 28311 if kind == "second_Q0" else 0
        for src in plan["sources"]:
            offsets[kind, src["index"]] = offset
            offset += src["dimension"]
        assert offset == (77639 if kind == "second_Q0" else 28311)
    assert [(s["index"],s["profile"],s["dimension"]) for s in plans["first"]["sources"]] == [
        (s["index"],s["profile"],s["dimension"]) for s in plans["second_Q1"]["sources"]]
    assert 256 * (PRIME - 1) ** 2 < 2**63

    def target_product(task):
        equation, target = task
        j = target["index"]
        rhs_parent = "first" if equation == 1 else "second_rhs"
        rhs_dir = output / rhs_parent / "rhs" / ("profile%02d" % j)
        stamp = read(rhs_dir / "result.json")
        inputs = read(rhs_dir / "manifest.json")
        assert stamp["target"] == j and inputs["target"] == target
        assert stamp["status"] == "complete_negative_Q%d_evaluations" % equation
        assert inputs["frame"]["prime"] == PRIME and inputs["frame"]["columns"] == stamp["columns"]
        rhs = vector(rhs_dir / "rhs.npy", stamp["sha256"], stamp["columns"])
        residual = -rhs
        standalone = -rhs.copy() if equation == 1 else None
        digests = []
        entries = 0
        for kind in (("first",) if equation == 1 else ("second_Q0", "second_Q1")):
            for src in plans[kind]["sources"]:
                if j not in src["targets"]:
                    continue
                i, n = src["index"], src["dimension"]
                folder = output / "compacted" / kind / "maps" / ("profile%02d" % i) / ("rows%05d_%05d" % (0,n))
                info = read(folder / "manifest.json")
                assert info["prime"] == PRIME and info["source"] == src
                assert info["start_row"] == 0 and info["stop_row"] == n
                assert info["targets"][str(j)] == inputs["frame"]
                base = "t%02d_%05d_%05d" % (j,0,n)
                receipt = read(folder / (base + ".json"))
                path = folder / (base + ".npy")
                assert (receipt["source"],receipt["target"],receipt["start"],receipt["stop"],receipt["columns"]) == (i,j,0,n,len(rhs))
                assert sha(path) == receipt["sha256"]
                matrix = np.load(path, allow_pickle=False, mmap_mode="r")
                assert matrix.shape == (n,len(rhs)) and matrix.dtype == np.dtype("<u4")
                start = offsets[kind,i]
                for lo in range(0,n,256):
                    hi = min(n,lo+256)
                    block = matrix[lo:hi].astype(np.int64)
                    assert np.all(block < PRIME)
                    residual = (residual + joint[start+lo:start+hi] @ block) % PRIME
                    if standalone is not None:
                        standalone = (standalone + first[start+lo:start+hi] @ block) % PRIME
                digests.append([kind,i,j,n,len(rhs),receipt["sha256"]])
                entries += matrix.size
        assert not residual.any()
        if standalone is not None:
            assert not standalone.any()
        return dict(equation=equation, target=j, coordinates=len(rhs), joint_residual_nonzero=0,
                    standalone_first_residual_nonzero=0 if standalone is not None else None,
                    matrix_entries=int(entries), rhs_sha256=stamp["sha256"], matrices=digests)

    tasks = [(eq,t) for eq,kind in ((1,"first"),(2,"second_Q0")) for t in plans[kind]["targets"]]
    with ThreadPoolExecutor(workers) as pool:
        records = list(pool.map(target_product, tasks))
    widths = [sum(r["coordinates"] for r in records if r["equation"] == eq) for eq in (1,2)]
    assert widths == [34901,65406]
    return dict(group="SU(3)", level=30, prime=PRIME, joint_coordinates=sum(widths),
                first_coordinates=widths[0], second_coordinates=widths[1], joint_residual_nonzero=0,
                standalone_first_residual_nonzero=0, matrix_entries=sum(r["matrix_entries"] for r in records),
                packed_matrices=sum(len(r["matrices"]) for r in records), targets=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="Existing production directory with runtime/ and output/")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    result = replay(args.root, args.workers)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k:v for k,v in result.items() if k != "targets"}))
