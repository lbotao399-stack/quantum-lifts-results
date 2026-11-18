"""Verify the published SU(3) solution coordinates, stage records and pivot choices.

The residual counts in the solve records were computed on the cluster. This
verifier needs NumPy and does not reconstruct the production map matrices.
"""
import hashlib
import gzip
import io
import json
import struct
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

import numpy as np


def verify(receipt_path):
    archive_bytes = Path(receipt_path).read_bytes()
    assert hashlib.sha256(archive_bytes).hexdigest() == "b4b4e1df18bfa3666d5c0f9b09443d06e0087a41a78fcb742078f8e378a9557e"
    with ZipFile(io.BytesIO(archive_bytes)) as archive:
        manifest = json.loads(archive.read("SHA256SUMS.json"))
        assert set(archive.namelist()) == set(manifest) | {"SHA256SUMS.json"}
        for name, digest in manifest.items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == digest
        first_prime = json.loads(archive.read("first_prime_stages.json"))
        original_pivots = []
        stage_hashes = []
        pivot_indices = []
        summaries = []
        solutions = []
        free_coordinates = []

        for folder, count in (("first_solve_compact", 13), ("second_solve_compact", 37)):
            result = json.loads(archive.read(folder + "/result.json"))
            stages = json.loads(archive.read(folder + "/stages.json"))
            assert result["prime"] == 1000033
            assert len(stages) == count
            assert [int(s["folder"].split("target")[1]) for s in stages] == result["processed_targets"]
            for stage in stages:
                pivots = stage["pivots"]
                assert stage["consistent"] and stage["rank"] == len(pivots)
                digest = hashlib.sha256(struct.pack(f"<{len(pivots)}q", *pivots)).hexdigest()
                original_pivots.append((stage["folder"], stage["rank"], digest))
                stage_hashes.append(stage["result_sha256"])
                pivot_indices.extend(pivots)
            assert stage_hashes == result["stage_certificates"]

            def array(name, sha_key):
                payload = archive.read(f"{folder}/{name}.npy")
                assert hashlib.sha256(payload).hexdigest() == result[sha_key]
                return np.load(io.BytesIO(payload), allow_pickle=False)

            particular = array("particular", "particular_sha256")
            free = array("free_variables", "free_variables_sha256")
            variables = result["variables"]
            assert len(set(pivot_indices)) == len(pivot_indices) == result["equation_rank"]
            assert all(0 <= i < variables for i in pivot_indices)
            assert np.array_equal(free, sorted(set(range(variables)) - set(pivot_indices)))
            assert len(free) == result["homogeneous_dimension"]
            assert particular.dtype.kind == "u" and particular.shape == (variables,)
            assert np.all(particular < result["prime"]) and np.all(particular[free] == 0)
            solutions.append(particular)
            free_coordinates.append(free)
            assert all(row["nonzero"] == 0 for row in result["residuals"])
            coordinates = Counter()
            profiles = set()
            for row in result["residuals"]:
                equation = row.get("equation", 1)
                profiles.add((equation, row["target"]))
                coordinates[equation] += row["coordinates"]
            expected = {(1, i) for i in range(13)}
            if count == 37:
                expected.update((2, i) for i in range(37))
            assert profiles == expected and len(profiles) == len(result["residuals"])
            summaries.append(dict(variables=variables, rank=len(pivot_indices),
                                  free_variables=len(free), coordinates_by_equation=dict(coordinates),
                                  nonzero_coordinates_reported=0,
                                  particular_sha256=result["particular_sha256"],
                                  free_variables_sha256=result["free_variables_sha256"]))
        assert all(row["consistent"] for row in first_prime)
        assert original_pivots == [(s["name"], s["rank"], s["pivots_sha256"]) for s in first_prime]
        assert summaries[0]["coordinates_by_equation"] == {1: 34901}
        assert summaries[1]["coordinates_by_equation"] == {1: 34901, 2: 65406}
        assert [(s["variables"], s["rank"], s["free_variables"]) for s in summaries] == [
            (28311, 18003, 10308), (77639, 50147, 27492)]
        first, joint = solutions
        changed = joint[:len(first)] != first
        first_free = free_coordinates[0]
        new_pivots = pivot_indices[summaries[0]["rank"]:]
        inherited_pivots = [i for i in new_pivots if i < len(first)]
        assert set(inherited_pivots) <= set(first_free)
        correction_comparison = dict(
            first_correction_coordinates_changed=int(np.count_nonzero(changed)),
            inherited_free_coordinates_changed=int(np.count_nonzero(changed[first_free])),
            new_pivots_in_first_correction=len(inherited_pivots),
            new_pivots_in_second_correction=len(new_pivots) - len(inherited_pivots),
        )
        return dict(matching_ordered_pivot_lists=len(original_pivots),
                    first_equation=summaries[0], coupled_equations=summaries[1],
                    correction_comparison=correction_comparison)


if __name__ == "__main__":
    data = Path(__file__).resolve().parents[1] / "data"
    path = data / "N4_SU3_j30_solution_receipts_20261118.zip"
    summary = verify(path)
    first = json.loads((data/"N4_SU3_j30_first_prime_residuals_20261118.json").read_text())
    second = json.loads((data/"N4_SU3_j30_independent_residuals_20261118.json").read_text())
    assert (first["prime"],second["prime"]) == (1000003,1000033)
    raw = gzip.decompress((data/"N4_SU3_j30_first_prime_solution.npy.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == first["solution_sha256"]
    coordinates = np.load(io.BytesIO(raw), allow_pickle=False)
    assert coordinates.dtype == np.uint32 and coordinates.shape == (77639,)
    assert np.all(coordinates < first["prime"])
    with ZipFile(path) as z:
        free = z.read("second_solve_compact/free_variables.npy")
        assert hashlib.sha256(free).hexdigest() == first["free_sha256"]
        assert not np.any(coordinates[np.load(io.BytesIO(free), allow_pickle=False)])
    left = {(r["order"],r["target"]): (r["coordinates"],r["nonzero"]) for r in first["residuals"]}
    right = {(r["equation"],r["target"]): (r["coordinates"],r["joint_residual_nonzero"]) for r in second["targets"]}
    assert left == right and len(left) == len(first["residuals"]) == len(second["targets"]) == 50
    assert all(v[1] == 0 for v in left.values())
    assert (first["variables"],first["rank"],first["free"]) == (77639,50147,27492)
    summary["first_prime_replay"] = dict(prime=first["prime"],solution_sha256=first["solution_sha256"],
        common_free_coordinates=first["free"],matching_target_receipts=len(left),
        first_coordinates=sum(n for (k,j),(n,v) in left.items() if k==1),
        second_coordinates=sum(n for (k,j),(n,v) in left.items() if k==2),nonzero_residuals=0)
    print(json.dumps(summary, indent=2))
