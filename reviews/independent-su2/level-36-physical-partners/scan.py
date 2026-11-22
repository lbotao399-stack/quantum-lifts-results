"""Compute a fixed physical charge in the full relative SU(2) Taylor complex."""
import argparse
import json
import os
from pathlib import Path
import subprocess


def run(engine, output, first_length, grades, prime, threads, seconds):
    env = dict(os.environ, PAR_THREADS=str(threads), A_VERB="0", A_CHECK="1")
    env.pop("A_NORHO", None)
    for grade in grades:
        weight = (6-grade,)*3 + (grade,)*2
        length = first_length-grade
        name = f"Y{first_length}-t{grade}-p{prime}"
        path = output/(name+".json")
        if not path.exists():
            seed = 850713 + 311*grade + 221*first_length + prime % 131
            args = [engine, "chain", "SU2", *weight, length, length, seed, prime]
            try:
                with (output/(name+".log")).open("w") as log:
                    proc = subprocess.run(list(map(str,args)), stdout=subprocess.PIPE,
                                          stderr=log, env=env, check=True, timeout=seconds)
            except subprocess.TimeoutExpired:
                print(name, "time limit", flush=True)
                continue
            data = json.loads(proc.stdout)
            assert (data["seed"], data["prime"]) == (seed, prime)
            assert data["gauge"] == "SU2" and data["sector"] == list(weight)
            row, = data["blocks"]
            assert row["Y"] == length and row["h"] == row["dim"]-row["rank_in"]-row["rank_out"]
            assert 0 <= row.get("rho",0) <= row["h"]
            path.write_text(json.dumps(data, indent=2)+"\n")
        print(name, json.loads(path.read_bytes())["blocks"][0], flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--first-length", type=int, default=9)
    p.add_argument("--prime", type=int, default=1000003)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--seconds", type=int, default=180)
    p.add_argument("--lease", action="store_true")
    p.add_argument("grades", type=int, nargs="+", choices=range(7))
    a = p.parse_args()
    if not __debug__ or not (a.threads > 0 and a.seconds > 0):
        p.error("enable assertions and use positive resource limits")
    a.output.mkdir(exist_ok=True, parents=True)
    bus = Path("~/.agent-bus/bin/agent-bus").expanduser()
    if a.lease:
        subprocess.run([bus,"lease",str(os.getpid()),"--agent","codex","--purpose",
                        "SU(2) level-36 physical partner ranks","--hours","2"],check=True)
    try:
        run(a.engine.resolve(),a.output,a.first_length,a.grades,a.prime,a.threads,a.seconds)
    finally:
        if a.lease:
            subprocess.run([bus,"release",str(os.getpid())],check=True)
