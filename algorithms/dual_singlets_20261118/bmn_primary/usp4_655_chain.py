"""Resolve the five odd lengths of the USp(4) level-32 BMN weight (6,5,5)."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def save(path, data):
    temporary = path.with_suffix(".part")
    temporary.write_text(json.dumps(data, indent=2)+"\n")
    temporary.replace(path)


def run(binary, output, threads, job):
    # Independent symplectic character dimensions in length order, 6 through 16.
    dimensions = (10, 84, 796, 3805, 9330, 13616, 13046, 8431, 3448, 739, 52)
    output.mkdir(parents=True, exist_ok=True)
    files = {}
    for length in (7,15,9,13,11):
        path = output/f"chain_6_5_5_{length}.json"
        seed = 761231+length
        prime = 2147483647
        if path.exists():
            record = json.loads(path.read_bytes())
        else:
            args = [str(binary.resolve()), "USp4", "6", "5", "5", str(length),
                    str(length), str(seed), str(threads), str(prime), "--probe-last"]
            record = json.loads(subprocess.check_output(args))
        assert record["group"] == "USp(4)" and record["weight"] == [6,5,5]
        assert record["prime"] == prime and record["seed"] == seed and record["terminal_probes_only"]
        block, = record["blocks"]
        assert (block["length"],block["dimension"]) == (length,dimensions[length-6])
        assert 0 <= block["incoming"] <= dimensions[length-7]
        assert 0 <= block["outgoing"] <= dimensions[length-5]
        assert block["incoming"]+block["outgoing"] <= block["dimension"]
        save(path, record)
        files[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        save(output/"receipt.partial.json", dict(job=job,files=files))
    (output/"receipt.partial.json").replace(output/"receipt.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--job", required=True)
    args = parser.parse_args()
    run(args.binary,args.output,args.threads,args.job)
