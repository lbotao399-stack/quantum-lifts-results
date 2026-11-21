"""Verify the Slurm receipts and combine the two finite-colour rank bounds."""
import gzip
import hashlib
import json
from math import prod
from pathlib import Path


def calculate(root):
    load = lambda p: json.loads(p.read_bytes())
    accounting = load(root / "jobs.json")
    for folder, name in [("first_action", "hairmem"), ("predecessors", "hairbound")]:
        source = root / folder
        receipt = load(source / "receipt.json")
        job, state, _, _, code = accounting[name].split("|")
        assert receipt["job"] == job and state == "COMPLETED" and code == "0:0"
        for file, digest in receipt["files"].items():
            assert hashlib.sha256((source / file).read_bytes()).hexdigest() == digest
    original = root.parents[1] / "so7_cartan_target_20261118/representatives.json.gz"
    rank = load(root / "first_action/first.json")
    assert hashlib.sha256(original.read_bytes()).hexdigest() == rank["source_sha256"]
    before = json.loads(gzip.decompress((root / "predecessors/previous.json.gz").read_bytes()))
    assert rank["weight"] == before["weight"] == [1, 1, 4, 3, 3]
    assert rank["prime"] == before["prime"] == 1000003
    assert (rank["level"], rank["source_terms"], rank["first_terms"], rank["first_denominator"]) == (26, 1668, 94926, 21)
    assert before["source_length"] == 8 and before["nine_letter_dimension"] == 25329
    panels = before["panels"]
    row_ids = [i for p in panels for i in p["rows"]]
    columns = [i for p in panels for i in p["columns"]]
    assert all(len(p["rows"]) == len(p["columns"]) for p in panels)
    assert len(set(row_ids)) == len(row_ids) == len(set(columns)) == len(columns) == before["closed_rank_lower_bound"]
    assert min(row_ids) >= 0 and max(row_ids) < before["formal_inputs"]
    assert min(columns) >= 0 and max(columns) < before["projected_coordinates"]
    det = prod(p["det"] for p in panels) % before["prime"]
    assert len(row_ids) == 12400 and det == 343790
    upper = before["nine_letter_dimension"] - len(row_ids)
    assert rank["classical_rank"] == upper == before["ten_letter_boundary_rank_upper_bound"]
    assert rank["augmented_rank"] == upper + 1 and not rank["first_correction_consistent_on_probes"]
    return dict(preceding_minor_order=len(row_ids), preceding_minor_residue=det,
                classical_rank_over_Q=upper, augmented_rank_over_Q=upper + 1,
                new_d1_pairings=1, new_fortuitous_generators=0)


if __name__ == "__main__":
    assert __debug__
    print(json.dumps(calculate(Path(__file__).resolve().parent), indent=2))
