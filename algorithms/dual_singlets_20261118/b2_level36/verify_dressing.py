"""Verify the rank bounds for the only possible lower long-root dressing."""
import gzip
import hashlib
import json
from math import prod
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from b2_schur_complex import multiplicity, transpose


def minor(source, job, length):
    job_id, state, _, _, exit_code = job.split("|")
    receipt = json.loads((source/"receipt.json").read_bytes())
    raw = (source/"image.json.gz").read_bytes()
    assert (state, exit_code, receipt["job"]) == ("COMPLETED", "0:0", job_id)
    assert receipt["files"] == {"image.json.gz": hashlib.sha256(raw).hexdigest()}
    record = json.loads(gzip.decompress(raw))
    assert (record["source_length"], record["group"], record["colours"]) == (length, "Sp", 4)
    assert record["weight"] == [0, 0, 8, 4, 4] and record["prime"] == 1000003
    rows, columns, determinants = [], [], []
    for panel in record["panels"]:
        assert len(panel["rows"]) == len(panel["columns"]) > 0
        assert 0 < panel["det"] < record["prime"]
        rows.extend(panel["rows"])
        columns.extend(panel["columns"])
        determinants.append(panel["det"])
    assert len(rows) == len(set(rows)) == len(set(columns)) == record["closed_rank_lower_bound"]
    assert len(rows) <= record["sampled_inputs"] <= record["formal_inputs"]
    assert min(rows) >= 0 and max(rows) < record["formal_inputs"]
    assert min(columns) >= 0 and max(columns) < record["projected_coordinates"]
    return dict(job=job_id, rank_lower=len(rows),
                minor_residue=prod(determinants) % record["prime"],
                target_dimension=record["target_dimension"])


def verify():
    folder = HERE/"dressing"
    jobs = json.loads((folder/"jobs.json").read_bytes())
    preceding, following = [minor(folder/str(y), jobs[str(y)], y) for y in (10, 11)]
    shape = transpose((8, 4, 4))
    dimensions = [multiplicity(shape, y) for y in (10, 11, 12)]
    index = json.loads(gzip.decompress((HERE/"index_columns.json.gz").read_bytes()))
    column, = [r for r in index if r["weight"] == [8, 4, 4]]
    assert dimensions == [column["dimensions"][str(y)] for y in (10, 11, 12)]
    assert dimensions[1:] == [preceding["target_dimension"], following["target_dimension"]]
    assert preceding["rank_lower"] + following["rank_lower"] == dimensions[1]
    lower = json.loads(gzip.decompress((HERE/"odd_through_j28.json.gz").read_bytes()))
    lower_blocks = {(tuple(w), y): (n, rin, rout) for w,y,n,rin,rout in lower["nonempty_blocks"]}
    wall = json.loads((HERE/"wall_restriction.json").read_bytes())
    assert wall["su2_odd_through_24"] == [[24, 7, [4,4,4], 1]]
    assert wall["long_wall_relative_to_mixed_coefficient"] == "-1/1728"
    candidates = []
    for spectator_degree in (2, 4, 6):
        active_weight, y = (10-spectator_degree, 4, 4), 13-spectator_degree
        if spectator_degree == 2:
            n, rin, rout = dimensions[1], preceding["rank_lower"], following["rank_lower"]
        else:
            n, rin, rout = lower_blocks[(active_weight, y)]
        assert n == rin+rout
        candidates.append(dict(spectator_degree=spectator_degree, weight=active_weight,
                               level=2*sum(active_weight), length=y, cochains=n,
                               incoming=rin, outgoing=rout, rational_cohomology=n-rin-rout))
    return dict(prime=1000003, adjacent_minors=[preceding, following],
                dressing_candidates=candidates, wall_relative_to_published_mixed_coefficient="-1/1728",
                published_input="Nonzero SU(4) mixed restriction, arXiv:2610.20340v2 Eq. (4.15)",
                independent_rational_USp4_level36_cycle=False)


if __name__ == "__main__":
    assert __debug__
    print(json.dumps(verify(), indent=2))
