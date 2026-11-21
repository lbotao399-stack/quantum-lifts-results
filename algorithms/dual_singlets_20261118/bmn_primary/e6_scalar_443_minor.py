"""Evaluate the odd split-Casimir scalar and the classical boundaries on the same E6 frames."""
import argparse
import json
from pathlib import Path

import numpy as np

from e6_casimir_443 import expansion, projections, split_casimir_values
from e6_scalar_boundaries import frame, m, trace_boundaries


HERE = Path(__file__).resolve().parent


def compute(trace_record, work):
    prime, seed, points = (trace_record[k] for k in ("prime", "seed", "points"))
    weight = (4, 4, 3)
    previous, = [row for row in trace_record["results"]
                 if row["weight"] == list(weight) and row["outer"] == -1]
    letters = [(0, 0, (i,)) for i in range(3)]
    samples = frame(letters, points, prime, seed)
    work.mkdir(parents=True, exist_ok=True)
    boundaries = trace_boundaries(weight, -1, samples, m.library(work, 27, prime), prime)
    boundary = m.minor(boundaries, prime)
    assert boundary == previous["boundary_minor"] and boundary["rank"] == 22

    triple = (0, 1, 2)
    polynomial = expansion(((0, 2), triple, triple, triple))
    projectors = projections(prime)
    inverse_two = pow(2, -1, prime)
    values = []
    for i in range(points):
        scalars = [samples.L[l][1][i, 0] for l in letters]
        twice, = split_casimir_values(scalars, [polynomial], prime, projectors, -1)
        values.append(twice * inverse_two % prime)
    augmented = np.vstack((boundaries, values))
    minor = m.minor(augmented, prime)
    assert minor["rank"] == previous["dimension"] == 23
    assert minor["rows"][-1] == len(boundaries)
    return dict(group="E6", weight=weight, level=22, length=11, outer_parity=-1,
                prime=prime, seed=seed, points=points, boundary_minor=boundary,
                scalar_values=values, augmented_minor=minor,
                augmented_matrix=augmented[np.ix_(minor["rows"], minor["columns"])].tolist())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", choices=(1000003, 1000033), type=int, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    record = json.loads((HERE/"scalar_443"/f"e6_trace_{args.prime}.json").read_bytes())
    result = compute(record, args.work)
    args.out.write_text(json.dumps(result, separators=(",", ":"))+"\n")
    print(result["prime"], result["augmented_minor"]["rank"], result["augmented_minor"]["determinant"])
