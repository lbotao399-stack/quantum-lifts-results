#!/usr/bin/env python3
"""Reproduce the two ordered stable derivatives requested in discussion D006.

Read fixed Git objects into a separate work directory; do not edit the source
checkout. Requires Python, git, Node.js, C++17, numpy, scipy, python-flint.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import base64
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import time

HERE = Path(__file__).resolve().parent
DOT = "e1d184d388c588c0a5a825a19ecbca4e24427242"
CODEX = "9a554a6dbf1792b4f3d552cd24bf1fc6aaaee6e9"
RUNTIME = "082be7070d207e2cdecfd0a91b03d5e014fd303d"
DUAL_PATH = "algorithms/dual_singlets_20261118"
RT_PATH = "algorithms/su3_independent_prime_20261118/runtime"


def git(repo, *args, **kwargs):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs).stdout


def materialize(repo, commit, prefix, target, select):
    try:
        git(repo, "cat-file", "-e", commit + "^{commit}")
    except subprocess.CalledProcessError:
        git(repo, "fetch", "--depth", "1", "--filter=blob:none", "origin", commit)
    rows = []
    for line in git(repo, "ls-tree", "-r", commit, prefix).decode().splitlines():
        header, name = line.split("\t", 1)
        mode, kind, sha = header.split()
        if kind == "blob" and select(name):
            rows.append((name, sha))
    assert rows
    git(repo, "fetch", "--no-tags", "--no-write-fetch-head", "origin", "--stdin",
        input=("\n".join(sorted({sha for name, sha in rows})) + "\n").encode())
    manifest = []
    for name, sha in rows:
        data = git(repo, "cat-file", "blob", sha)
        path = target / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        manifest.append({"commit": commit, "path": name, "git_blob": sha,
                         "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=HERE / "so11_replay_work")
    args = parser.parse_args()
    repo, work = args.repo.resolve(), args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    sources = work / "pinned_sources"
    manifest = materialize(repo, CODEX, DUAL_PATH, sources, lambda name:
        len(Path(name).parts) == 3 or
        (len(Path(name).parts) == 4 and "/orthogonal_glueballs/" in name) or
        "/orthogonal_glueballs/so11_j36/" in name)
    manifest += materialize(repo, RUNTIME, RT_PATH, sources, lambda name: True)
    manifest += materialize(repo, DOT, "reviews/dot-orthogonal-d3", sources,
                            lambda name: name.endswith(".md") or name.endswith("so11-j36-replay.tar.gz.b64"))
    encoded = sources / "reviews/dot-orthogonal-d3/so11-j36-replay.tar.gz.b64"
    archive = base64.b64decode(encoded.read_bytes())
    assert hashlib.sha256(archive).hexdigest() == "d55e71942d7e3f72db765b9f97c31b0e234b0afcd3b70b2e1b8fbf011587bc2b"
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tf:
        tf.extractall(work / "dot", filter="data")
    dot = work / "dot/so11-j36"
    for line in (dot / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert hashlib.sha256((dot / name.lstrip("*")).read_bytes()).hexdigest() == digest

    fresh = work / "fresh_banks"
    fresh.mkdir(exist_ok=True)
    executable = fresh / "bmngram"
    subprocess.run(["c++", "-O3", "-std=c++17", str(dot / "bmngram.cpp"), "-o", str(executable)], check=True)
    circuits = []
    for label, script in [("cycle", "seed_cycle.py"), ("path3", "seed_3path.py"), ("path3-alt", "seed_3path_alt.py")]:
        seed = subprocess.check_output([sys.executable, str(dot / script)], text=True).splitlines()
        assert seed[0] == "3 0"
        seed[0] = "6 0"
        seed = "\n".join(seed) + "\n"
        raw = subprocess.check_output([str(executable)], input=seed, text=True)
        (fresh / (label + ".raw")).write_text(raw)
        circuits.append({"seed": label, "input": seed, "rows": len(raw.splitlines()),
                         "matching_counter": raw.splitlines()[-1],
                         "sha256": hashlib.sha256(raw.encode()).hexdigest()})
    probe = HERE / "so11_derivative_probes.js"
    archived = json.loads(subprocess.check_output(["node", str(probe), str(dot)], text=True))
    regenerated = json.loads(subprocess.check_output(["node", str(probe), str(dot), str(fresh)], text=True))
    assert archived == regenerated
    independent = json.loads(subprocess.check_output([
        sys.executable, str(HERE / "so11_derivative_review.py"),
        "--archived", str(dot), "--fresh", str(fresh),
        "--output", str(work / "so11_derivative_review.json")], text=True))
    assert independent["passed"]

    os.environ["N4_TRACE_RUNTIME"] = str(sources / RT_PATH / "work/sp3-native.nosync")
    so11 = sources / DUAL_PATH / "orthogonal_glueballs/so11_j36"
    sys.path.insert(0, str(so11))
    sys.path.insert(0, str(so11.parent))
    import terminal_primitive as terminal
    import align_pair
    terminal_data = terminal.check_terminal(True)
    alignment = align_pair.main(sources / "reviews/dot-orthogonal-d3")
    rows = []
    for row in regenerated["output"]:
        codex = terminal.leading_derivative(row["query"], terminal.PfaffianCoefficient(3, 0))
        dot_value = Fraction(row["signed_dot_derivative"])
        rows.append({"ordered_traces": row["query"], "dot": str(dot_value), "codex": str(codex),
                     "closed_combination_derivative": str(9 * dot_value + 25600 * codex)})
    determinant = Fraction(rows[0]["dot"]) * Fraction(rows[1]["codex"]) - Fraction(rows[1]["dot"]) * Fraction(rows[0]["codex"])
    assert determinant == 2654208000
    result = {"passed": True, "sources": {"dot": DOT, "codex": CODEX, "runtime": RUNTIME},
              "ordered_derivative_rows": rows, "formal_defect_rank": 2,
              "determinant": str(determinant), "fresh_endpoints_match_archived_banks": True,
              "circuits": circuits, "dot_projection_and_returns": regenerated,
              "independent_fraction_check": independent,
              "alignment": alignment, "terminal": terminal_data,
              "scope": "Chosen stable derivatives only. No claim that the two cohomology classes are independent, or that the closed combination is nontrivial modulo the full Q0(I).",
              "seconds": time.monotonic() - started}
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    (work / "INPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("terminal", "circuits", "dot_projection_and_returns")}, indent=2))


if __name__ == "__main__":
    main()
