"""Test the shared level-32 correction on a noncommuting three-colour evaluation."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

from su3_level32_walls import DATA, algebra, block_first, block_reduce
from replay_third_action import basis_values
from verify import add_scaled


def check(data, certificate):
    info = json.loads((data / "order1.json").read_bytes())
    correction = {block_reduce.to_fort(w):Fraction(c) for w,c in certificate["correction"]}
    boundary = algebra.differential(correction)
    rows = defaultdict(dict)
    direct = defaultdict(dict)
    digest = hashlib.sha256()
    with gzip.open(data / "Q1.jsonl.gz","rb") as stream:
        for line in stream:
            digest.update(line)
            value, raw = json.loads(line)
            word = tuple(bytes(map(block_first.code,tr)) for tr in raw)
            t = bytes(sorted(b"".join(word)))
            direct[t][word] = Fraction(value)
    assert digest.hexdigest() == info["output_sha256"] and sum(map(len,direct.values())) == info["direct_terms"]
    for full,c in boundary.items():
        word = tuple(bytes(tr) for tr in block_reduce.from_fort(full))
        rows[bytes(sorted(b"".join(word)))][word] = c
    tested = 0
    for t in sorted(set(direct) | set(rows), key=lambda t:(len(direct[t])+len(rows[t]),t)):
        terms = direct[t].copy()
        add_scaled(terms, rows[t], 1)
        if not terms:
            continue
        target = tuple(sorted(map(block_first.letter,t)))
        def evaluation(poly):
            return basis_values({target:{block_reduce.to_fort(w):c for w,c in poly.items()}}, 346712)
        tested += 1
        probe = evaluation(terms)
        value = Fraction(probe["contents"][0]["coefficient"])
        if value:
            first = evaluation(direct[t])["contents"][0]["coefficient"]
            classical = evaluation(rows[t])["contents"][0]["coefficient"]
            assert value == Fraction(first)+Fraction(classical)
            return dict(first_image_sha256=digest.hexdigest(), correction_terms=len(correction),
                        formal_boundary_terms=len(boundary), nonempty_contents_tested=tested,
                        Q1_source=first, Q0_correction=classical, residual=str(value), evaluation=probe)
    raise ValueError("no residual found at this evaluation")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, default=DATA)
    p.add_argument("--certificate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    result=check(a.data,json.loads(a.certificate.read_bytes()))
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result["Q1_source"],result["Q0_correction"],result["residual"])
