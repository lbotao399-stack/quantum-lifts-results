"""Verify native minors in the alternative physical polarization of the rank-three Levi diagrams."""

import argparse
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess

from levi_source_minors import PRIME_SEEDS, validate


def verify(engine=None, case=None, threads=1):
    base = Path(__file__).resolve().parent
    archive = json.loads(
        gzip.decompress((base / "levi_alternate_minors.json.gz").read_bytes())
    )
    table = json.loads((base / "levi_polarizations.json").read_bytes())
    schur = json.loads((base / "levi_relation_spaces.json").read_bytes())
    entries = archive["entries"]
    assert set(entries) == set(archive["receipts"])
    if case and case not in entries:
        raise ValueError(f"unknown case: {case}")
    rows = []
    for name, native in entries.items():
        receipt = archive["receipts"][name]
        job = archive["jobs"][receipt["job"]]
        assert job["state"] == "COMPLETED" and job["exit_code"] == "0:0"
        assert (
            sha256(
                (json.dumps(native, separators=(",", ":")) + "\n").encode()
            ).hexdigest()
            == receipt["sha256"]
        )
        group = native["gauge"]
        weight = native["sector"]
        prime = native["prime"]
        mapping = {(4, 4, 0, 2, 2): (4, 4, 2, 2), (4, 4, 4, 1, 1): (4, 4, 4, 2)}
        shape = mapping[tuple(weight)]
        diagram = next(
            d for d in table[group]["diagrams"] if tuple(d["shape"]) == shape
        )
        assert diagram["selected"] == weight and diagram["branching_multiplicity"] == 1
        check = next(
            d["independent_physical_control"]
            for d in schur[group]
            if tuple(d["shape"]) == shape
        )
        assert check["weight"] == weight
        dims = [native[k] for k in ("parent", "cochains", "next")]
        assert dims == diagram["selected_sizes"] == check["reconstructed"]
        ranks = validate(native, group, shape, prime, weight)
        if engine is not None and (case is None or case == name):
            raw = subprocess.check_output(
                [
                    str(engine),
                    *map(str, weight),
                    str(native["length"]),
                    str(prime),
                    str(PRIME_SEEDS[prime]),
                    group,
                ],
                env=dict(os.environ, PAR_THREADS=str(threads)),
            )
            assert json.loads(raw) == native
        rows.append(
            dict(
                job=job["slurm_id"],
                group=group,
                shape=shape,
                weight=weight,
                prime=prime,
                adjacent_dimensions=dims,
                ranks=ranks,
                ordinary_cohomology_upper_bound=ranks[0] - ranks[1] - ranks[2],
                determinants=[
                    native[k]["det"]
                    for k in ("incoming", "outgoing", "outgoing_and_Cartan")
                ],
            )
        )
    assert set(archive["jobs"]) == {r["job"] for r in rows}
    return sorted(rows, key=lambda r: (r["shape"], r["group"], r["prime"]))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", type=Path)
    p.add_argument("--case")
    p.add_argument("--threads", type=int, default=1)
    a = p.parse_args()
    if not __debug__ or a.threads < 1 or (a.case and not a.engine):
        p.error("enable assertions and use positive threads; --case requires --engine")
    print(
        json.dumps(
            verify(a.engine.resolve() if a.engine else None, a.case, a.threads),
            indent=2,
        )
    )
