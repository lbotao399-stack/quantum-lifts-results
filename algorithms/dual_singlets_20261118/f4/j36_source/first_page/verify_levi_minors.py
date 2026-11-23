"""Verify the regular-scalar subgroup records and optionally regenerate their native minors."""

import argparse
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess

from levi_source_minors import PRIME_SEEDS, branches, taylor_weight, validate


def check(engine=None, case=None, threads=1):
    archive = json.loads(gzip.decompress(Path(__file__).with_name("levi_minors.json.gz").read_bytes()))
    entries = archive["entries"]
    summary = {}
    expected = set()
    if case and case not in entries:
        raise ValueError(f"absent certificate: {case}")
    for group in ("SU3", "Sp2", "SO7", "Sp3"):
        shapes = branches(2) if group in ("SU3", "Sp2") else [
            shape for shape in branches(1) if len(shape) <= 3 or shape == (4, 4, 3, 1)]
        rows = []
        for shape in shapes:
            witnesses = []
            for prime, seed in PRIME_SEEDS.items():
                name = "-".join(map(str, (group, *shape, prime))) + ".json"
                expected.add(name)
                entry = entries[name]
                encoded = (json.dumps(entry, separators=(",", ":")) + "\n").encode()
                assert sha256(encoded).hexdigest() == archive["scheduler_sha256"][name]
                ranks = validate(entry, group, shape, prime)
                assert ranks[-1] == 0
                witnesses.append(ranks)
                if engine is not None and (case is None or case == name):
                    args = [str(engine), *map(str, taylor_weight(shape)), str(sum(shape) - 1),
                            str(prime), str(seed), group]
                    fresh = subprocess.check_output(args, env=dict(os.environ, PAR_THREADS=str(threads)))
                    assert json.loads(fresh) == entry
            assert witnesses[0] == witnesses[1]
            rows.append(dict(shape=shape, ranks=witnesses[0]))
        summary[group] = rows
    assert set(entries) == set(archive["scheduler_sha256"]) == expected
    return dict(jobs=archive["jobs"], two_prime_blocks=sum(map(len, summary.values())), groups=summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path)
    p.add_argument("--case", help="regenerate only this JSON member")
    p.add_argument("--threads", type=int, default=1)
    a = p.parse_args()
    if not __debug__ or a.threads < 1 or (a.case and not a.engine):
        p.error("enable assertions, supply positive threads, and pass an engine when selecting a case")
    result = check(a.engine.resolve() if a.engine else None, a.case, a.threads)
    print(json.dumps(result, indent=2))
