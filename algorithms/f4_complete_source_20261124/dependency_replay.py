"""Replay the explicitly scoped public prerequisites of the sixth/seventh calculations."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
from types import SimpleNamespace

from survival_inputs import verify


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("outgoing", "product-ward", "seventh-first-page", "seventh-intermediate", "seventh-colour"))
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = verify(args.inputs)
    started = perf_counter()
    commit = "c778cd7c3ceb72ad8082d791b2aee5f0df082164"
    root = args.inputs / commit
    source = root / "algorithms/dual_singlets_20261118/f4/j36_source"
    if args.stage == "outgoing":
        commit = "3b30b8ff242977debd90581c12f4a767d7875c6c"
        root = args.inputs / commit
        directory = root / "reviews/independent-f4/sixth-glueball-classical-rank"
        sys.path.insert(0, str(directory))
        from verify_classical import verify as classical
        result = classical(SimpleNamespace(samples=[], full=False))
        scope = "Fresh rational characters, explicit Cartan Gram complement and determinants of both saved complete 703-square matrices. No historical native-column rebuilder or colour-subspace map is invoked."
        result["native_column_positions_rebuilt_at_each_prime"] = []
    elif args.stage == "product-ward":
        directory = source / "first_page"
        sys.path.insert(0, str(directory))
        from glueball_factor_ward import calculate
        result = calculate(0)
        expected = json.loads((directory / "glueball_factor_ward.json").read_text())
        assert json.loads(json.dumps(result)) == expected
        result = {"full_record_matches_pinned": True, "terminal_dimension": 55,
                  "rational_terminal_cokernel": result["upper_bound"],
                  "rational_reduced_determinant": result["harmonic_control"]["determinant"],
                  "rational_unreduced_determinant": result["joint"]["determinant"]}
        scope = "Fresh native rational reconstruction of all 25 extremal coordinates and the complete 8-dimensional boundary test needed by the sixth-power classical-product theorem."
    elif args.stage == "seventh-first-page":
        directory = source / "eight_albert_review"
        sys.path.insert(0, str(directory))
        from first_page_tower import Descent, check_power
        packet = json.loads((directory / "tower_probes.json").read_text())
        result = {"universal_scalar_check": Descent().check_psi(),
                  "powers": [check_power(7, packet["7"])]}
        assert result["powers"][0]["target_pairing"] == -194893931232000
        scope = "Fresh native rational first-page certificate for all seventh-power incoming profiles: 2704 scalar pairs, 1378 symmetric pairs, 1117 differentiated-scalar entries and 18948 mixed-ghost components. This does not decide the second differential."
    else:
        directory = source / "focused_powers"
        sys.path.insert(0, str(directory))
        if args.stage == "seventh-intermediate":
            from seventh_intermediate import calculate
            result = calculate()
            expected = json.loads((directory / "seventh-intermediate.json").read_text())
            assert json.loads(json.dumps(result)) == expected
            result["homogeneous_first_corrections"] = 23895 - result["intermediate_classical_rank"]
            scope = "Fresh native 68-square homogeneous Ward matrix and complete terminal character. This proves rational rank 2959 of the full first-correction map; it does not decide the seventh-power second map."
        else:
            from seventh_scalar_ward import run
            result = run()
            expected = json.loads((directory / "seventh-scalar-ward.json").read_text())
            assert json.loads(json.dumps(result)) == expected
            scope = "Exact character and stored scalar-weight rank consistency replay, with fresh rational invariant determinant. The stored 54 scalar-weight ranks are not all regenerated. No complete second-map decision is made."
    result["scope"] = scope
    result["public_commit"] = commit
    result["stage"] = args.stage
    result["elapsed_seconds"] = perf_counter() - started
    result["input_manifest_sha256"] = sha256(Path(__file__).with_name("survival_input_manifest.json").read_bytes()).hexdigest()
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("scalar_quotient",)}), flush=True)


if __name__ == "__main__":
    main()
