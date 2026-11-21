"""Bound the finite-colour classical image by evaluating exact formal cycles."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import mixed_incoming as m
from dual_charge_bands import coordinates
from mixed33_exact import apply, exact
from orthosymplectic_classical_page import minor
from trace_descent import highest_weights


def run(previous_path, target_path, out, workers):
    out.mkdir(parents=True, exist_ok=True)
    previous, target = [json.loads(gzip.decompress(p.read_bytes()))
                        for p in (previous_path, target_path)]
    weight = (3, 3, 5, 0, 0)
    for page, length in ((previous, 8), (target, 9)):
        assert (tuple(page["weight"]), page["length"], page["group"], page["colours"],
                page["prime"]) == (weight, length, "Sp", 6, m.P)
    source = highest_weights(weight, 8)
    assert set(source) == {3} and set(previous["centre_frames"]) == {"3"}
    rows = source[3]
    images = [m.q(row) for row in rows]
    # Exact elimination uses the formal multi-trace words as independent columns.
    independent, _, relations = exact(images)
    assert all(not apply(relation, images) for relation in relations)
    cycles = [apply(relation, rows) for relation in relations]
    alphabet = sorted({letter for content in m.contents(weight, 8) for letter in content})
    probes = previous["centre_frames"]["3"]["samples"]
    library = m.library(out, 6, m.P)
    evaluated = coordinates(cycles, 3, alphabet, probes, 6, "Sp", library, workers)
    witness = minor(evaluated, m.P)
    upper = min(len(independent), previous["cochains"] - witness["rank"])
    lower = target["incoming"]["rank"]
    assert upper >= lower == previous["outgoing"]["rank"]
    record = dict(weight=weight, prime=m.P, previous_sha256=hashlib.sha256(previous_path.read_bytes()).hexdigest(),
                  target_sha256=hashlib.sha256(target_path.read_bytes()).hexdigest(),
                  formal_domain=len(rows), formal_image_rank=len(independent),
                  formal_kernel_dimension=len(cycles), finite_domain=previous["cochains"],
                  finite_kernel_minor=witness, finite_image_lower_bound=lower,
                  finite_image_upper_bound=upper,
                  independent_formal_cycles=[sorted(relations[i].items()) for i in witness["rows"]])
    m.write(out/"boundary_bounds.json.gz", record)
    print("formal image and kernel", len(independent), len(cycles), "finite kernel minor",
          witness["rank"], "image interval", lower, upper, flush=True)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("previous", "target", "out"):
        parser.add_argument("--"+option, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    assert __debug__
    run(args.previous, args.target, args.out, args.workers)
