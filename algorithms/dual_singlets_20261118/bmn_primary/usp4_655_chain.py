"""Compute both adjacent classical ranks at each odd length of a rank-two BMN column."""
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
    (7,5,3): (83,893,4273,10637,14783,12131,6037,1792,248),
    (7,5,4): (52,986,5732,16884,29157,31504,21748,9356,2256,232),
    (8,5,2): (99,1063,3723,6713,7045,4299,1421,181),
}

SEEDS = {(6,5,5):761231, (6,5,4):793611, (7,4,4):829137,
         (7,5,3):863451, (7,5,4):889373, (8,5,2):919363}


def run(binary, output, threads, job, weight=(6,5,5), group="USp4"):
    if group == "USp4":
        dimensions = dict(enumerate(CHARACTERS[weight], max(weight)))
        seed_base = SEEDS[weight]
        group_argument = ["USp4"]
        lengths = [y for y in (7,15,9,13,11) if y in dimensions]
    else:
        # The independent G2 colour character has zero cochains at length five.
        if group != "G2" or weight != (5,5,5):
            raise ValueError("No independent character for this gauge group and flavour weight")
        dimensions = dict(enumerate((10,82,393,1387,3030,3834,2843,1300,416,85), 6))
        seed_base = 941771
        group_argument = []
        lengths = sorted((y for y in dimensions if y % 2), key=lambda y: dimensions[y])
    output.mkdir(parents=True, exist_ok=True)
    files = {}
    assert set(lengths) == {y for y in dimensions if y % 2}
    for length in lengths:
        path = output/("chain_"+"_".join(map(str,(*weight,length)))+".json")
        seed = seed_base+length
        prime = 2147483647
        if path.exists():
            record = json.loads(path.read_bytes())
        else:
            args = [str(binary.resolve()), *group_argument, *map(str,weight), str(length),
                    str(length), str(seed), str(threads), str(prime), "--probe-last"]
            record = json.loads(subprocess.check_output(args))
        assert record.get("group", "G2") == ("USp(4)" if group == "USp4" else group)
        assert record["weight"] == list(weight)
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
    parser.add_argument("--group", choices=("USp4", "G2"), default="USp4")
    args = parser.parse_args()
    choices = CHARACTERS if args.group == "USp4" else {(5,5,5)}
    if tuple(args.weight) not in choices or not __debug__ or args.threads < 1:
        parser.error("Choose a weight with an independent character table and enable assertions")
    run(args.binary,args.output,args.threads,args.job,tuple(args.weight),args.group)
