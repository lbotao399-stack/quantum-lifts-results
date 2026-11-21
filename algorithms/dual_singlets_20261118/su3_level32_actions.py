"""Compute the two direct supercharge images of the level-32 SU(3) BMN source."""
import argparse
from fractions import Fraction
import gzip
import hashlib
from itertools import islice
import json
from math import lcm
from pathlib import Path
import sys

from su3_level32_replay import DATA, LEVEL32_SHA
from second_action import cartan, read_source, save_rows
from verify import add_scaled


def calculate(runtime, source_path, order, out):
    sys.path.insert(0, str(runtime.resolve()))
    import su_quantum_integer_v1 as quantum
    import su_trace_adapter_v1 as algebra

    source = read_source(source_path, source_sha=LEVEL32_SHA, terms=8640, sector=(5,5,16,8,8))
    entries = iter(sorted(source.items()))
    image = {}
    controls = []
    for start in range(0, len(source), 128):
        chunk = dict(islice(entries, 128))
        if not controls:
            for i, (word, c) in enumerate(chunk.items(), start):
                value = quantum.reference({word:c}, order, 3)
                if value:
                    assert quantum.act({word:c}, order, 3) == value
                    controls.append(dict(index=i, terms=len(value)))
                    break
        add_scaled(image, quantum.act(chunk, order, 3), 1)
    assert controls and image and algebra.conjugate(image) == image
    assert {algebra.fw.sector_of(w) for w in image} == {(4,4,16,8,8)}
    assert {sum(map(len,w)) for w in image} == {12-order}
    diagonal = cartan(image)
    out.mkdir(parents=True, exist_ok=True)
    rows = ([str(image[w]), w] for w in sorted(image))
    digest = save_rows(out/f"Q{order}.jsonl.gz", rows)
    record = dict(group="SU(3)", level=32, order=order, source_sha256=LEVEL32_SHA,
                  source_terms=len(source), direct_terms=len(image), output_length=12-order,
                  doubled_output_charges=[4,4,16,8,8], output_sha256=digest, input_controls=controls,
                  common_denominator=lcm(*(Fraction(c).denominator for c in image.values())),
                  diagonal_trace_terms=diagonal["formal_diagonal_trace_terms"],
                  Cartan_polynomial=diagonal["terms"], diagonal_alphabet=diagonal["alphabet"])
    with gzip.open(out/f"Q{order}.jsonl.gz","rb") as f:
        assert hashlib.file_digest(f, "sha256").hexdigest() == digest
    (out/f"order{order}.json").write_text(json.dumps(record, indent=2)+"\n")
    print(order, len(image), diagonal["formal_diagonal_trace_terms"], len(diagonal["terms"]), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=DATA/"O32_C_even.json.gz")
    p.add_argument("--runtime", type=Path, default=DATA.parents[2]/"su3_independent_prime_20261118/runtime/work/sp3-native.nosync")
    p.add_argument("--order", type=int, choices=(1,2), required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    calculate(args.runtime, args.source, args.order, args.out)
