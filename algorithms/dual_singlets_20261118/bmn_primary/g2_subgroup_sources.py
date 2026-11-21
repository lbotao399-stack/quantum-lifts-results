"""Restrict the two possible odd G2 generator weights to its compact rank-one and rank-two subgroups."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
DATA = HERE / "subgroup_sources"
sys.path.insert(0, str(HERE.parent))
from bmn_primary_index import cochains
from generator_bounds import dressings
from rank_two_bmn_product_replay import highest, witness
from su3_bmn_parity_products import labelled
from trace_descent import dimensions
from verify import receipt


PUBLIC = {
    "su2/bmn_hw_cohomology_j44.json": "332f24c8c86a750005c60121592df804855ffb25",
    "su3/hw_cohomology_j32_records/su3_c_6_6_4.jsonl": "6ada104198a0cf995a90ff675f6683c8e2c7e563",
}
WEIGHTS = ((6,6,4), (6,6,6))


def upstream_records(directory):
    contents = {}
    for name, digest in PUBLIC.items():
        raw = (directory / name).read_bytes()
        assert hashlib.sha1(f"blob {len(raw)}\0".encode()+raw).hexdigest() == digest
        contents[name] = raw
    su2, su3 = contents.values()
    factors = []
    for column in json.loads(su2)["columns"]:
        weight = column["lam"]
        if weight[0] > 6:
            continue
        assert column["j"] == 2*sum(weight) and column["all_checks"]
        for row in column["sectors"]:
            assert row["check"] and 0 <= row["h"] <= row["m"]
            if row["m"]:
                factors.append([weight, sum(weight)-row["twoJ"], row["m"], row["h"]])
    unitary = {}
    for row in map(json.loads, su3.splitlines()):
        assert row["lam"] == [6,6,4] and row["N"] == 3 and row["check"]
        length = 16-row["n"]
        assert length not in unitary and 0 <= row["r_up"] <= row["m"]
        unitary[length] = [length, row["m"], row["r_up"]]
    return dict(commit="4256f207986029540060fb44dd494746828f65ce", blobs=PUBLIC,
                su2=sorted(factors), su3=[unitary[i] for i in sorted(unitary)])


def direct_sum_bounds(factors):
    levels = cochains(18,"SU(2)")[1]
    dims = {(w,y):n for y,rows in levels.items() for w,(n,) in rows.items() if w[0] <= 6}
    upper = {((0,0,0),0):1}
    for weight, length, dim, h in factors:
        key = tuple(weight), length
        assert key not in upper and dims[key] == dim and 0 <= h <= dim
        upper[key] = h
    assert dims.keys() == upper.keys()
    rows = []
    for target in WEIGHTS:
        products = Counter()
        chains = Counter()
        for (w,i), d in dims.items():
            for v, n in dressings(target,w).items():
                for j in range(max(v), sum(v)+1):
                    if (v,j) in dims:
                        chains[i+j] += n*d*dims[v,j]
                        products[i+j] += n*upper[w,i]*upper[v,j]
        for y in range(max(target)|1, sum(target)+1, 2):
            rows.append(dict(weight=target, length=y, two_factor_cochains=chains[y],
                             two_factor_cohomology_upper=products[y], principal_cochains=dims.get((target,y),0),
                             principal_cohomology_upper=upper.get((target,y),0)))
    assert {(r["weight"],r["length"]):r["two_factor_cohomology_upper"] for r in rows
            if r["two_factor_cohomology_upper"]} == {((6,6,4),9):2, ((6,6,6),9):6}
    return rows


def long_root_sector(rows):
    weight = (6,6,4)
    dims = cochains(16,"SU(3)")[1]
    incoming = {y:rank for y,_,rank in rows}
    assert {y for y,v in dims.items() if weight in v} == set(incoming)
    for y, dim, rank in rows:
        assert dim == sum(dims[y][weight]) and rank + incoming.get(y+1,0) <= dim
    upper = {y:sum(dims[y][weight])-incoming[y]-incoming.get(y+1,0) for y in (11,13)}
    assert upper == {11:1, 13:0}
    product = json.loads((DATA/"su3_664_product11.json").read_bytes())
    factors, = product["inputs"]
    factors = tuple(tuple(spec) for spec in factors)
    assert sum(sum(spec[:5]) for spec in factors) == product["length"] == 11
    polynomial = highest(labelled(factors), weight)
    minor = witness([polynomial])
    assert json.loads(json.dumps(minor)) == product["minor"]
    assert (minor["denominators"], minor["signed_integer_determinant"]) == ([5],"18")
    return [dict(length=y, total_cohomology_upper=upper[y], incoming=incoming[y], outgoing=incoming[y+1],
                 cochains=sum(dims[y][weight]), odd_colour_product_rank=int(y==11), even_cohomology_upper=0)
            for y in (11,13)]


def even_sector(path, weight):
    record = json.loads(path.read_bytes())
    assert record["sector"] == [*weight,0,0] and record["gauge"] == "SU3" and record["mode"] == "even_BMN_minors"
    length = record["length"]
    table = DATA/("su3_even_"+"".join(map(str,weight))+".txt")
    expected = {tuple(map(int,line.split())) for line in table.read_text().splitlines()}
    for y, field in ((length-1,"parent"),(length,"cochains"),(length+1,"target_cochains")):
        confirmed = dimensions((0,0,*weight), y, 3, "SU+")
        assert {(y,f,n) for f,n in confirmed.items()} == {k for k in expected if k[0] == y}
        assert sum(confirmed.values()) == record[field]
    for field, source_dim in (("incoming",record["parent"]),("outgoing",record["cochains"])):
        minor = record[field]
        assert 0 < minor["det"] < record["prime"] and minor["rank"] <= source_dim
        assert len(set(minor["columns"])) == len(set(minor["rows"])) == minor["rank"]
        assert max(minor["columns"]) < source_dim
    return dict(weight=weight, length=length, cochains=record["cochains"], incoming=record["incoming"]["rank"],
                outgoing=record["outgoing"]["rank"], prime=record["prime"],
                determinants=[record[field]["det"] for field in ("incoming","outgoing")])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, help="data directory at the pinned Fortuity commit")
    parser.add_argument("--extract", type=Path, help="save the compact upstream rank table")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__ or (args.extract and not args.upstream):
        parser.error("Enable assertions and supply the pinned data directory for extraction")
    pinned = DATA/"upstream_ranks.json"
    if args.upstream:
        records = upstream_records(args.upstream)
        if args.extract:
            args.extract.write_text(json.dumps(records,separators=(",",":"))+"\n")
        elif records != json.loads(pinned.read_bytes()):
            raise ValueError("Extracted rank table differs from pinned records")
    else:
        records = json.loads(pinned.read_bytes())
    assert records["blobs"] == PUBLIC
    result = dict(direct_sum=direct_sum_bounds(records["su2"]), long_root_664=long_root_sector(records["su3"]))
    receipt("7658790")
    result["long_root_666"] = even_sector(HERE/"7658790/even_15.json", (6,6,6))
    controls = [even_sector(DATA/name,(4,4,4)) for name in ("su3_even_444_control.json","su3_even_444_positive.json")]
    independent = json.loads((DATA/"su3_444_independent.json").read_bytes())
    assert independent["length"] == 9 and independent["sector"] == [4,4,4,0,0]
    assert (independent["even_in"],independent["even_out"]) == (controls[0]["incoming"],controls[0]["outgoing"])
    retained = [result["long_root_666"],*controls]
    assert [r["cochains"]-r["incoming"]-r["outgoing"] for r in retained] == [0,0,1]
    result["controls"] = controls
    args.out.write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
