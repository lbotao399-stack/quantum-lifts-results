"""Compute the odd lengths of a selected USp(4) BMN column."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def save(path, data):
    temporary = path.with_suffix(".part")
    temporary.write_text(json.dumps(data, indent=2)+"\n")
    temporary.replace(path)


CHARACTERS = {
    (6,5,5): (10,84,796,3805,9330,13616,13046,8431,3448,739,52),
    (6,5,4): (15,226,1622,5637,10889,12711,9285,4187,1077,116),
    (7,4,4): (38,758,3254,6962,9053,7432,3638,916,75),
}


def run(binary, output, threads, job, weight=(6,5,5)):
    # Dimensions in ascending letter length, beginning at the largest flavour count.
    dimensions = dict(enumerate(CHARACTERS[weight],max(weight)))
    output.mkdir(parents=True, exist_ok=True)
    files = {}
    for length in (7,15,9,13,11):
        path = output/("chain_"+"_".join(map(str,(*weight,length)))+".json")
        seed = {(6,5,5):761231,(6,5,4):793611,(7,4,4):829137}[weight]+length
        prime = 2147483647
        if path.exists():
            record = json.loads(path.read_bytes())
        else:
            args = [str(binary.resolve()), "USp4", *map(str,weight), str(length),
                    str(length), str(seed), str(threads), str(prime), "--probe-last"]
            record = json.loads(subprocess.check_output(args))
        assert record["group"] == "USp(4)" and record["weight"] == list(weight)
        assert record["prime"] == prime and record["seed"] == seed and record["terminal_probes_only"]
        block, = record["blocks"]
        assert (block["length"],block["dimension"]) == (length,dimensions[length])
        assert 0 <= block["incoming"] <= dimensions.get(length-1,0)
        assert 0 <= block["outgoing"] <= dimensions.get(length+1,0)
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
    parser.add_argument("--weight",type=int,nargs=3,default=(6,5,5))
    args = parser.parse_args()
    if tuple(args.weight) not in CHARACTERS or not __debug__:
        parser.error("Choose a weight with an independent character table and enable assertions")
    run(args.binary,args.output,args.threads,args.job,tuple(args.weight))
