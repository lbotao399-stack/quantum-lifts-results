"""Authenticate the scoped packet and replay its complete small rational identities.

The public order-3017 minor is an attributed lemma. This command checks its
recorded arithmetic, not its native entries. The literal 3017-by-50 extension
and the full bilinear theorem have separate verification scopes below.
"""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
from flint import fmpz_mat

from survival_inputs import pinned_files, verify

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text())


def replay(script, inputs, stored):
    with TemporaryDirectory(prefix="f4-packet-") as directory:
        result_path = Path(directory) / "result.json"
        subprocess.run([sys.executable, "-B", str(HERE / script),
                        "--inputs", str(inputs), "--out", str(result_path)],
                       check=True, capture_output=True, text=True)
        result = json.loads(result_path.read_text())
    assert result == read(stored), script
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = verify(args.inputs)
    _, records = pinned_files()
    public = args.inputs / "c778cd7c3ceb72ad8082d791b2aee5f0df082164"
    focused = public / "algorithms/dual_singlets_20261118/f4/j36_source/focused_powers"
    product = replay("verify_closed_product.py", focused, "closed_product_rational.json")
    bilinears = replay("verify_closed_bilinears.py", args.inputs, "closed_bilinears_rational.json")
    tree = fmpz_mat(bilinears["twice_classical_matrix"])
    kernel, nullity = tree.nullspace()
    assert tree.rank() == 20 and nullity == 40 and (tree * kernel).is_zero()
    assert (fmpz_mat([bilinears["first_quantum_row"]]) * kernel).is_zero()

    extension = read("literal_albert_extension_1009.json")
    compressed = (HERE / "literal_albert_extension_1009.u16.gz").read_bytes()
    raw = gzip.decompress(compressed)
    assert sha256(compressed).hexdigest() == extension["gzip_sha256"]
    assert sha256(raw).hexdigest() == extension["residue_u16le_sha256"]
    assert len(raw) == 2 * 3017 * 50 and extension["shape"] == [3017, 50]
    assert extension["complete"] and extension["completed_rows"] == 3017
    assert len(extension["points"]) == len(set(extension["points"])) == 3017
    bank = np.frombuffer(raw, dtype="<u2").reshape(3017, 50)
    assert int(bank.max()) < extension["prime"] == 1009
    for name, expected in extension["implementation_sha256"].items():
        assert sha256((HERE / name).read_bytes()).hexdigest() == expected
    controls = read("albert_integer_controls_1009.json")
    for control in controls["checks"]:
        position = extension["points"].index(control["point"])
        assert bank[position].tolist() == control["residues"]
        assert control["verified_entries"] == 50 and control["classical_split_nonzeros"] == 0

    historical_minors = {}
    selection = None
    for prime, old, complement, determinant in ((1009, 833, 737, 449), (1013, 615, 346, 60)):
        certificate = json.loads((focused / f"sixth-product-{prime}.json").read_text())
        assert (certificate["previous_minor"], certificate["complementary_minor"], certificate["determinant"]) == (old, complement, determinant)
        assert old * complement % prime == determinant != 0
        assert certificate["rank"] == 3017 and certificate["previous_rank"] == 3002
        identity = certificate["new_points"], certificate["new_columns"]
        if selection is None:
            selection = identity
        assert identity == selection and certificate["new_columns"][-1] == 49
        historical_minors[str(prime)] = determinant

    outgoing = read("outgoing_rank_replay.json")
    assert outgoing["classical_rank"] == 703 and outgoing["classical_kernel_dimension"] == 3017
    assert outgoing["native_column_positions_rebuilt_at_each_prime"] == []
    assert outgoing["minor_determinants"] == {"1009": 1005, "1013": 734}
    assert outgoing["cartan_value"] == 1728 * 5184 == 8957952
    assert outgoing["dimensions"] == [10007, 3720, 704]
    assert 3016 + 1 + 703 == 3720
    result = {
        "schema": "f4-sixth-scoped-packet-verification-v1",
        "field": "Q, with modular nonzero minors supplying rational lower bounds",
        "verified_public_files": len(records),
        "numerical_input_files": len(manifest["files"]),
        "input_manifest_sha256": sha256((HERE / "survival_input_manifest.json").read_bytes()).hexdigest(),
        "proof_dependencies_sha256": sha256((HERE / "proof_dependencies.json").read_bytes()).hexdigest(),
        "fresh_rational_product": {"columns": 96, "equations": 172, "rank": product["rational_rank"],
                                   "nonzero_coefficients": 87, "closed_spectator_pairs": 45},
        "fresh_complete_bilinear_identity": {"columns": 60, "classical_rank": 20,
                                             "complete_kernel_dimension": 40, "quantum_rank_on_kernel": 0},
        "authenticated_literal_extension_entries": int(bank.size),
        "literal_extension_prime": 1009,
        "attributed_complete_kernel_minors": historical_minors,
        "fresh_native_complete_kernel_minor_rebuild": False,
        "outgoing_saved_matrix_determinants": outgoing["minor_determinants"],
        "complete_kernel_dimension": 3720 - 703,
        "classical_boundaries_plus_closed_product": [3016, 1],
        "all_homogeneous_first_corrections": 4109 - 633,
        "complete_second_map_rank": 0,
        "decision_basis": "Audited published full-kernel rank theorem, exact direct closed-factor proof, matched quantum-product chain identity and the published full first-page Ward identity.",
        "fresh_verification_scope": "Complete small rational product and bilinear identities; authenticated full literal extension; saved outgoing determinants. The historical 3002-column block and full 3017 determinant are not freshly rebuilt.",
        "C04_complete_predecessor_basis": False,
        "D04_complete_seventh_second_map": False,
        "physical_instanton_protection_inferred": False,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
