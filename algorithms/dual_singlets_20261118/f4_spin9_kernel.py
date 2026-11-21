"""Compute the classical rank bounds on the full-rank Spin(9) subgroup of F4."""
import argparse
import hashlib
import json
import os
from itertools import permutations
from pathlib import Path
import shutil
import subprocess

from charge_maps import _multiplicity, content_patterns
from su3_levi_kernel import check, seed_for
from trace_descent import dimensions


HERE = Path(__file__).resolve().parent
HOOKS = HERE / "taylor_hooks"


def compile_kernel(work):
    directory = work / "source"
    directory.mkdir(exist_ok=True)
    pins = json.loads((HOOKS / "native/sha256.json").read_bytes())
    for name, digest in pins.items():
        source = HOOKS / "native" / name
        if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected native source {name}")
        shutil.copyfile(source, directory / name)
    for name in ("unitary_parity.cpp", "su2_cartan_kernel.cpp"):
        shutil.copyfile(HOOKS / name, directory / name)
    patch = HERE / "f4/spin9_kernel.patch"
    subprocess.run(["patch", "-s", "-p1", "-i", str(patch)], cwd=directory, check=True)
    executable = work / "spin9-kernel"
    subprocess.run(["c++", "-O3", "-std=c++17", "-pthread", str(directory / "su2_cartan_kernel.cpp"),
                    "-o", str(executable)], check=True)
    return executable


def character_table(width, path):
    patterns = set()
    for length in (3*width-2, 3*width-1, 3*width):
        for order in permutations(range(3)):
            weight = (0, 0, *(width + order[i] - i for i in range(3)))
            for bosons, fermions, _ in content_patterns(weight, length):
                entries = [(count, 0) for count in bosons]
                entries.extend((count, int(count > 1)) for count in fermions)
                patterns.add(tuple(sorted(entries)))
    with path.open("w") as output:
        print("SO9", file=output)
        for pattern in sorted(patterns):
            bosons = tuple(count for count, odd in pattern if not odd)
            fermions = tuple(count for count, odd in pattern if odd)
            value = _multiplicity(bosons, fermions, 0, 9)
            print(len(pattern), *(v for entry in pattern for v in entry), value, file=output)
    return len(patterns)


def validate(record, width, prime):
    weight = (width, width, width, 0, 0)
    length = 3*width-1
    cochains, upper, kernel = check(record, weight, length, prime, "SO9")
    sizes = [sum(dimensions((0, 0, width, width, width), y, 9, "SO").values())
             for y in (length-1, length, length+1)]
    if sizes != [record["parent"], cochains, record["next"]]:
        raise ValueError("Frame dimensions disagree with the orthogonal character")
    return dict(weight=list(weight), length=length, dimensions=sizes, prime=prime,
                incoming_at_least=record["incoming"]["rank"], outgoing_at_least=record["outgoing"]["rank"],
                outgoing_and_Cartan_at_least=record["outgoing_and_Cartan"]["rank"],
                classical_cohomology_at_most=upper, Cartan_kernel_at_most=kernel)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, choices=(4, 5), required=True)
    parser.add_argument("--prime", type=int, choices=(1000003, 1000033, 2147483647), default=2147483647)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--records", type=Path, help="Replay the saved rank and character bounds")
    args = parser.parse_args()
    if not __debug__ or args.threads < 1:
        parser.error("The replay requires assertions and a positive thread count")
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.records:
        record = json.loads(args.records.read_bytes())
    else:
        executable = compile_kernel(work)
        table = work / "character.txt"
        n = character_table(args.width, table)
        weight = (args.width, args.width, args.width, 0, 0)
        length = 3*args.width - 1
        seed = seed_for(weight, length)
        print(f"Spin(9), length {length}, {n} character patterns", flush=True)
        command = [str(executable), *map(str, (*weight, length, args.prime, seed)), "SO9"]
        env = dict(os.environ, N4_B4_PATTERNS=str(table), PAR_THREADS=str(args.threads))
        record = json.loads(subprocess.check_output(command, env=env))
        (work / "minors.json").write_text(json.dumps(record, separators=(",", ":")) + "\n")
    result = validate(record, args.width, args.prime)
    (work / "bounds.json").write_text(json.dumps(result, indent=2) + "\n")
    print(result, flush=True)
