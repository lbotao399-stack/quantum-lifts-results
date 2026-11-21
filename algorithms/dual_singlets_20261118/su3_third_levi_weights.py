"""Bound the level-36 third target on the two-colour Levi in three one-flavour components."""
import argparse
import gzip
import json
from pathlib import Path

from su3_levi_branching import branches
from su3_levi_kernel import compute, gzip_json, key, resolve


PRIME = 2147483647


def cases():
    result = {}
    # Removing a vertical strip from (6,6,1,1,1) leaves one of these two-row diagrams.
    for flavour, plus, minus in ((3, 6, 6), (4, 6, 5), (5, 5, 5)):
        weight = (flavour, 0, 0, plus, minus)
        for length, parity in ((8, -1), (9, 1), (10, -1)):
            result[f"{flavour}_{length}"] = dict(weight=weight, length=length, colour_parity=parity,
                                                  active=branches(weight, length, parity))
    return result


def read_records(paths):
    records = {}
    for path in paths:
        batch = json.loads(gzip.decompress(path.read_bytes()))
        for label, record in batch.items():
            if label in records and records[label] != record:
                raise ValueError(f"incompatible records for {label}")
            records[label] = record
    return records


def missing(plan, recorded):
    return {label: dict(active=[row for row in case["active"]
                                if key(row["weight"], row["length"]) not in recorded])
            for label, case in plan.items()}


def bounds(plan, records):
    sectors = resolve(plan, records, PRIME)["sectors"]
    result = {}
    for length in (8, 9, 10):
        terms = [sectors[f"{flavour}_{length}"] for flavour in (3, 4, 5)]
        result[str(length + 2)] = dict(target_length=length, colour_parity=terms[0]["parity"],
                                       classical_upper=min(t["classical_dimension_upper"] for t in terms),
                                       cartan_kernel_upper=min(t["cartan_kernel_upper"] for t in terms))
    return dict(prime=PRIME, partition=[6, 6, 1, 1, 1], sector_bounds=sectors,
                bounds_by_source_length=result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--records", type=Path, action="append", required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--engine", type=Path)
    p.add_argument("--raw", type=Path)
    p.add_argument("--workers", type=int, default=2)
    args = p.parse_args()
    if not __debug__ or bool(args.engine) != bool(args.raw) or args.workers < 1:
        p.error("use --raw with --engine and enable assertions")
    plan = cases()
    known = read_records(args.records)
    if args.engine:
        fresh = compute(missing(plan, known), args.engine, args.raw, args.workers, PRIME)
        gzip_json(args.out, fresh)
    else:
        args.out.write_text(json.dumps(bounds(plan, known), indent=2) + "\n")
