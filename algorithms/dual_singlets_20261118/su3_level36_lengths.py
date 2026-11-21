"""Bound individual rational cohomology dimensions in the level-36 SU(3) BMN column."""
import argparse
import hashlib
import json
from pathlib import Path

from su3_bmn_parity_replay import replay
from trace_descent import dimensions


HERE = Path(__file__).resolve().parent / "taylor_hooks/su3_bmn_sources"


def read(path):
    return json.loads(path.read_bytes())


def certified_files(folder, accounting):
    job, state, _, _, status = accounting.split("|")
    receipt = read(folder/"receipt.json")
    if receipt["job"] != job or state != "COMPLETED" or status != "0:0":
        raise ValueError(f"incomplete job for {folder.name}")
    for name, digest in receipt["files"].items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"file differs from its A8 receipt: {name}")
    return [folder/name for name in receipt["files"]]


def calculate(folder):
    weight = (8, 5, 5)
    prime = 2147483647
    products = replay(folder/"weight_8_5_5.json")
    assert json.loads(json.dumps(products)) == read(folder/"replay_weight_8_5_5.json")
    jobs = read(folder/"jobs.json")
    certified_files(folder/"descendants", jobs["descendants"])
    descendants = read(folder/"descendants/descendants.json")
    even_descendant = read(folder/"descendants/level34.json")["C_even_length11_raised"]
    cochains = {y: [sum(dimensions((0,0,*weight), y, 3, c).values()) for c in ("SU+", "SU-")]
                for y in range(max(weight), sum(weight)+1)}
    lower = {y: [0, 0] for y in cochains}
    for y, data in products["ordinary_diagonal_minors"].items():
        lower[y][y % 2] = data["rank"]
    for y, key in [(10, "C_odd_length10_triple_raising"), (12, "C_odd_length12_bilinear_product")]:
        assert descendants[key]["rational_class_nonzero"] and int(descendants[key]["pairing"])
        lower[y][1] += 1
    assert even_descendant["rational_class_nonzero"] and int(even_descendant["pairing"])
    lower[11][0] = 1
    if "image_11" in jobs:
        descendants32 = read(folder.parent/"su3_level32/descendant.json")
        assert descendants32["descendant_weight"] == list(weight) and descendants32["descendant_length"] == 11
        (a,b),(c,d) = [[int(value) for value in row] for row in descendants32["pairings"]]
        determinant = a*d-b*c
        assert b == 0 and d == int(even_descendant["pairing"])
        assert determinant == int(descendants32["determinant"]) == int(read(folder/"quartic.json")["determinant"])
        assert determinant
        lower[11][0] = 2
    image = {min(cochains)-1: [0, 0], max(cochains): [0, 0]}
    sources = [("classical_images", "image")]
    if "image_11" in jobs:
        sources.append(("length11_image", "image_11"))
    for directory, job in sources:
        for path in certified_files(folder/directory, jobs[job]):
            record = read(path)
            y = record["length"]
            if (record["mode"], record["gauge"], record["sector"], record["prime"]) != (
                    "colour_image", "SU3", [*weight, 0, 0], prime):
                raise ValueError(f"incorrect image at {path}")
            assert y not in image and record["cochains"] == sum(cochains[y])
            image[y] = []
            for i, name in enumerate(("even", "odd")):
                minor = record[name]
                rank = minor["rank"]
                assert 0 <= rank <= min(cochains[y][i], cochains[y+1][i])
                assert 0 < minor["det"] < prime
                for key, n in (("rows", record["target_evaluations"]), ("columns", record["cochains"])):
                    selected = minor[key]
                    assert len(selected) == len(set(selected)) == rank and all(0 <= j < n for j in selected)
                image[y].append(rank)
    result = []
    for y, total in cochains.items():
        if y-1 not in image or y not in image:
            continue
        upper = [d-i-o for d,i,o in zip(total, image[y-1], image[y])]
        assert all(a <= b for a,b in zip(lower[y], upper))
        result.append(dict(length=y, cochains=total, incoming_lower=image[y-1], outgoing_lower=image[y],
                           rational_lower=lower[y], rational_upper=upper))
    return dict(group="SU(3)", level=36, weight=weight, parities=["even", "odd"], prime=prime, lengths=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--folder", type=Path, default=HERE)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    data = calculate(args.folder)
    args.out.write_text(json.dumps(data, indent=2) + "\n")
