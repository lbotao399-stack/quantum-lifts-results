"""Compute the rational right-hand side of the first level-34 correction equation."""
import argparse
import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import islice
from pathlib import Path

from second_action import BMN_LETTERS, HERE, SOURCE_SHA, cartan, read_source, save_rows


def calculate(runtime, source_path, output):
    sys.path.insert(0, str(runtime.resolve()))
    import su_trace_adapter_v1 as algebra
    import su_quantum_integer_v1 as quantum
    from sp_full_profiles import profile
    from verify import add_scaled

    source = read_source(source_path)
    examples = {}
    rhs = {}
    entries = iter(sorted(source.items()))
    while chunk := dict(islice(entries, 128)):
        for word, c in chunk.items():
            key = any(BMN_LETTERS[-1] in trace for trace in word)
            if key not in examples:
                test = quantum.act({word: c}, 1, 3)
                if test:
                    assert test == quantum.reference({word: c}, 1, 3)
                    examples[key] = len(test)
        add_scaled(rhs, quantum.act(chunk, 1, 3), -1)
    assert set(examples) == {False, True} and rhs and algebra.conjugate(rhs) == rhs
    assert {algebra.fw.sector_of(w) for w in rhs} == {(5, 5, 16, 10, 8)}
    assert {sum(map(len, w)) for w in rhs} == {11}
    by_profile = Counter(profile(tuple(l for t in w for l in t)) for w in rhs)
    diagonal = cartan(rhs)
    output.mkdir(parents=True, exist_ok=True)
    sha = save_rows(output / "negative_Q1.jsonl.gz", ([str(rhs[w]), w] for w in sorted(rhs)))
    record = dict(source_sha256=SOURCE_SHA, source_terms=len(source), sign=-1, order=1,
                  formal_trace_terms=len(rhs), rhs_uncompressed_sha256=sha,
                  profiles=[dict(profile=p, terms=n) for p, n in sorted(by_profile.items())],
                  common_denominators=sorted({Fraction(v).denominator for v in rhs.values()}),
                  Cartan_formal_terms=diagonal["formal_diagonal_trace_terms"],
                  Cartan_polynomial_terms=len(diagonal["terms"]),
                  direct_single_word_checks=[examples[False], examples[True]])
    (output / "first_action.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k != "profiles"}, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--runtime", type=Path, default=HERE.parent / "su3_independent_prime_20261118/runtime/work/sp3-native.nosync")
    p.add_argument("--source", type=Path, default=HERE / "O34_C_even.json.gz")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    calculate(args.runtime, args.source, args.output)
