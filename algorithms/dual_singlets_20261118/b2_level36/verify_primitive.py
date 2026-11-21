"""Locate the odd cohomology in the (10,4,4) column and check its lower BMN sources."""
import gzip
import hashlib
import json
import sys
from pathlib import Path

from primitive_sources import candidates
from verify_dressing import verify as verify_dressing

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from b2_schur_complex import multiplicity, transpose


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def job_files(folder, accounting):
    job, state, _, _, exit_code = accounting.split("|")
    receipt = read(folder/"receipt.json")
    if (job, state, exit_code) != (receipt["job"], "COMPLETED", "0:0"):
        raise ValueError(f"incomplete job at {folder.name}")
    for name, expected in receipt["files"].items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"job output does not match its receipt: {name}")
    return receipt["files"]


def length_bounds():
    weight = (10, 4, 4)
    data = read(HERE/"rational_products.json")
    shape = transpose(weight)
    dims = {int(k): sum(s.values()) for k, s in data["cochains_by_length_and_fermion_degree"].items()}
    assert tuple(data["weight"]) == weight
    assert dims == {y: multiplicity(shape, y) for y in range(10, 19)}
    even = sum(n for y, n in dims.items() if y % 2 == 0)
    odd = sum(dims.values()) - even
    products = {int(y): m["rank"] for y, m in data["even_product_minors"].items()
                if int(m["signed_integer_determinant"]) != 0}
    assert (even, odd) == (data["even_cochains"], data["odd_cochains"])
    ranks = {9: 0, 18: 0}
    for y in (10, 11):
        image = read(HERE/f"length{y}.json.gz")
        minor = image["rank_minor"]
        assert image["source_dimension"] == dims[y] and 0 < minor["determinant"] < image["prime"]
        assert len(set(minor["rows"])) == len(set(minor["columns"])) == minor["rank"]
        ranks[y] = minor["rank"]
    path = HERE/"length12.json.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == read(HERE/"length12_receipt.json")["files"]["image.json.gz"]
    image = read(path)
    assert image["source_length"] == 12 and image["target_dimension"] == dims[13]
    rows, columns = [], []
    for panel in image["panels"]:
        assert len(panel["rows"]) == len(panel["columns"]) and 0 < panel["det"] < image["prime"]
        rows.extend(panel["rows"])
        columns.extend(panel["columns"])
    assert len(rows) == len(set(rows)) == len(set(columns)) == image["closed_rank_lower_bound"]
    ranks[12] = len(rows)
    folder = HERE/"high_lengths"
    assert set(job_files(folder, read(folder/"jobs.json")["chain"])) == {"chain.json"}
    chain = read(folder/"chain.json")
    assert chain["sector"] == [*weight, 0, 0] and chain["gauge"] == "Sp2"
    assert [b["Y"] for b in chain["blocks"]] == [15, 16, 17]
    for b in chain["blocks"]:
        y = b["Y"]
        assert dims[y] == b["dim"] == b["rank_in"] + b["rank_out"]
        assert ranks.get(y-1, b["rank_in"]) == b["rank_in"]
        ranks[y-1], ranks[y] = b["rank_in"], b["rank_out"]
    fixed = {}
    for y in sorted(dims.keys() - {13, 14}):
        upper = dims[y] - ranks[y-1] - ranks[y]
        assert upper == products.get(y, 0)
        fixed[y] = upper
    offset13 = dims[13] - ranks[12]
    offset14 = dims[14] - ranks[14]
    assert offset13 == offset14
    forced = sum(products.values()) - (even - odd)
    assert forced == products[14] == data["odd_cohomology_lower_bound"]
    return dict(weight=weight, level=36, cochains=dims, fixed_rational_cohomology=fixed,
                length13_and_length14_common_offset_before_subtracting_rank13=offset13,
                length13_cohomology_lower=forced, length14_cohomology_lower=forced,
                length13_image_upper=offset13-forced)


def check():
    column = length_bounds()
    dressing = verify_dressing()["dressing_candidates"][0]
    measured = {(tuple(dressing["weight"]), dressing["length"]): dressing["rational_cohomology"]}
    folder = HERE/"parents"
    for name, accounting in read(folder/"jobs.json").items():
        for filename in job_files(folder/name, accounting):
            row = read(folder/name/filename)
            weight, y = tuple(row["weight"]), row["length"]
            assert row["group"] == "USp(4)" and row["level"] == 2*sum(weight)
            dimension = multiplicity(transpose(weight), y)
            upper = dimension-row["incoming"]-row["outgoing"]
            assert dimension == row["cochains"] and upper == row["rational_cohomology_upper"] >= 0
            assert (weight, y) not in measured
            measured[weight, y] = upper
    channels = candidates()
    for channel in channels:
        channel["rational_parent_upper"] = measured.get((tuple(channel["weight"]), channel["length"]))
    products = [r for r in channels if "product" in r["operation"]]
    raising = [r for r in channels if r not in products]
    assert len(products) == 7 and all(r["rational_parent_upper"] == 0 for r in products)
    result = dict(column=column, lower_channels=channels, product_subspace_dimension=0,
                  indecomposable_modulo_products_lower=column["length13_cohomology_lower"])
    if all(r["rational_parent_upper"] is not None for r in raising):
        bound = sum(r["multiplicity"]*r["rational_parent_upper"] for r in raising)
        result.update(raising_image_upper=bound,
                      modulo_products_and_BMN_raising_lower=max(0, result["indecomposable_modulo_products_lower"]-bound))
    return result


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("verification requires enabled assertions")
    print(json.dumps(check(), indent=2))
