"""Check the exact permitted public packet; --write records its final hashes."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from survival_inputs import pinned_files

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "publication_manifest.json"
ALLOWED = (
    "D03_PROOF.md", "D04_STATUS.md", "METHOD_AUDIT.md", "PUBLICATION.md", "README.md",
    "albert_integer_controls_1009.json", "albert_tensor_runner.cpp", "albert_tensor_runner.py",
    "character_backend/b2_schur_complex.py", "character_backend/f4_singlet_characters.py",
    "character_inputs.json", "character_inventory.py", "closed_bilinears_rational.json",
    "closed_product_rational.json", "dependency_replay.py", "first_correction_inventory.json",
    "literal_albert_extension_1009.json", "literal_albert_extension_1009.u16.gz",
    "outgoing_rank_replay.json", "packet_verification.json", "product_ward_replay.json",
    "proof_dependencies.json", "publication.py", "rebuild_albert_extension.py", "requirements.txt",
    "second_source_inventory.json", "seventh_first_page_replay.json", "seventh_intermediate_replay.json",
    "seventh_scalar_ward_replay.json", "seventh_ward_complete_inventory.json", "seventh_ward_inventory.py",
    "stopped_route_provenance.json", "survival_input_manifest.json", "survival_inputs.py",
    "verify_closed_bilinears.py", "verify_closed_product.py", "verify_packet.py",
)


def current():
    actual = {str(path.relative_to(HERE)) for path in HERE.rglob("*")
              if path.is_file() and "__pycache__" not in path.parts and path != MANIFEST}
    if actual != set(ALLOWED):
        raise ValueError({"missing_allowed_files": sorted(set(ALLOWED)-actual),
                          "unlisted_files": sorted(actual-set(ALLOWED))})
    files = []
    for relative in sorted(ALLOWED):
        raw = (HERE / relative).read_bytes()
        files.append({"path": relative, "bytes": len(raw), "sha256": sha256(raw).hexdigest()})
    _, inputs = pinned_files()
    inputs = [{key: item[key] for key in ("commit", "path", "bytes", "sha256")}
              for item in sorted(inputs, key=lambda item:(item["commit"], item["path"]))]
    return {
        "schema": "f4-permitted-public-packet-v1", "date": "2026-11-24",
        "manifest_self_digest_excluded": True,
        "files": files, "public_inputs": inputs,
        "public_input_repository": "lbotao399-stack/quantum-lifts-results",
        "input_file_count": len(inputs), "input_bytes": sum(item["bytes"] for item in inputs),
        "result_scope": "Complete sixth-power second-map zero from the audited public full-rank theorem and the fresh exact direct closed-factor proof.",
        "fresh_native_complete_minor_rebuilt": False,
        "literal_extension_entries": 150850,
        "stopped_route_source_included": False,
        "C04_complete_predecessor_basis": False,
        "D04_complete_seventh_second_map": False,
        "physical_instanton_matrix_element": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = current()
    if args.write:
        MANIFEST.write_text(json.dumps(result, indent=2) + "\n")
    else:
        if json.loads(MANIFEST.read_text()) != result:
            raise ValueError("public packet differs from its final manifest")
    print(json.dumps({"manifest_sha256": sha256(MANIFEST.read_bytes()).hexdigest(),
                      "permitted_packet_files": len(result["files"]),
                      "public_input_files": result["input_file_count"],
                      "public_input_bytes": result["input_bytes"],
                      "fresh_native_complete_minor_rebuilt": False}))


if __name__ == "__main__":
    main()
