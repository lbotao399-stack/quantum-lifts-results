"""Multiply the five rational corrections and compare the two archived extractions."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import json
from pathlib import Path

import block_first
import block_reduce
from su3_levi_coefficient_solve import WEIGHTS, formal, polynomial
from verify import add_scaled


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def decode(rows):
    return {tuple(bytes.fromhex(t) for t in traces): Fraction(c) for traces, c in rows}


def compare_previous(target, path):
    old = {tuple(map(tuple, traces)): c for c, traces in map(json.loads, gzip.open(path, "rb"))}
    converted = Counter()
    for word, value in target.items():
        rewritten = tuple(tuple(block_first.code(formal.ALPH[k]) for k in t) for t in word)
        sign, representative = block_reduce.canonical(rewritten)
        if sign:
            converted[representative] -= value * sign
    assert {w: c for w, c in converted.items() if c} == old
    return len(old)


def replay(folder):
    data = read(folder / "coefficients.json.gz")
    assert data["negative_input_sha256"] == formal.INPUT and data["denominator"] == formal.DEN
    results = {}
    for name, weight in WEIGHTS.items():
        target = decode(data["coefficients"][name])
        record = read(folder / (name + ".json.gz"))
        assert record["channel"] == name and record["weight"] == list(weight)
        correction = decode(record["correction_after_spectator"])
        sign = -1 if name in ("three_phi1_minus", "three_phi1_plus") else 1
        residual = {w: c / formal.DEN for w, c in target.items()}
        for word, c in correction.items():
            add_scaled(residual, formal.q(word), sign*c)
        actual = polynomial(target)
        assert len(actual) == record["nonzero_input_monomials"] and not polynomial(residual)
        result = dict(input_terms=len(target), input_monomials=len(actual), correction_terms=len(correction),
                      nonzero_residual_monomials=0)
        if name in ("three_phi1_minus", "three_phi1_phi3"):
            old = "lambda_minus" if name == "three_phi1_minus" else "phi3"
            root = Path(block_first.__file__).resolve().parent / "block_first"
            result["matching_previous_terms_after_reversal"] = compare_previous(target, root / (old + ".jsonl.gz"))
        results[name] = result
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    args.out.write_text(json.dumps(replay(args.folder), indent=2) + "\n")
