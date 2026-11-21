"""Solve the level-32 first correction on the three scalar collision loci."""
import argparse
import json
from pathlib import Path
import sys

from su3_level32_replay import DATA


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("extract", "solve", "verify"))
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--scalar", choices=(1, 2, 4), type=int)
    parser.add_argument("--power", type=int)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    runtime = DATA.parents[2] / "su3_independent_prime_20261118/runtime/work/sp3-native.nosync"
    sys.path.insert(0, str(runtime))
    import sp_environment
    from collision_first import extract
    from collision_solve import solve
    from replay_collision import run

    record = json.loads((args.data / "order1.json").read_bytes())
    assert (record["order"], record["source_terms"], record["output_length"]) == (1, 8640, 11)
    options = dict(rhs_sha=record["output_sha256"], denominator=record["common_denominator"])
    if args.mode == "extract":
        extract(args.data / "Q1.jsonl.gz", args.work, input_words=record["direct_terms"],
                sign_of_input=-1, **options)
        return
    manifest = json.loads((args.work / "extraction.json").read_bytes())
    assert manifest["sign_applied_to_input"] == -1
    if args.mode == "verify":
        result = run(args.work, args.work, **options)
        (args.work / "verified.json").write_text(json.dumps(result, indent=2) + "\n")
        return
    if args.scalar is None or args.power is None:
        parser.error("solve requires --scalar and --power")
    solve(args.work, args.scalar, args.power, args.work / f"case-{args.scalar}-{args.power}.json",
          mu=(1, 1, 7, 3, 3), **options)


if __name__ == "__main__":
    main()
