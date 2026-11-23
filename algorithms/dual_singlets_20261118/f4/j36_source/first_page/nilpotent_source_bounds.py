"""Replay all transverse scalar degrees through seven for the selected F4 source."""

import argparse
import gzip
import json
from pathlib import Path

from flint import fmpq, fmpq_mat
from nilpotent_arcs import calculate as straight
from nilpotent_cubic import kernel
from nilpotent_curved_arcs import calculate as curved
from nilpotent_quadratic import array, calculate as quadratic, pivots
from nilpotent_stabilizer import reflection_word, source_action
from nilpotent_tangent import calculate as tangent, simple_coordinates


def invariant_bound(q, short, order, first, second):
    if not second["global_source_symbol_upper_bound"]:
        return dict(bound=0)
    _, fr, mons, inc, out = quadratic(q, short, order, True)
    n = len(mons)
    equations = out.tolist()
    for row in first["independent_evaluations"] + second["additional_evaluations"]:
        coeff = dict(row["evaluation"])
        equations.append([fmpq(coeff.get(i, 0)) for i in range(n)])
    vectors = kernel(fmpq_mat(equations))
    _, keep = pivots(inc)
    boundaries = array(n, [[inc[i, j] for i in range(n)] for j in keep])
    assert (
        vectors.ncols() - boundaries.ncols()
        == second["global_source_symbol_upper_bound"]
    )
    identity = fmpq_mat([[int(i == j) for j in range(n)] for i in range(n)])
    native = []
    for w in reflection_word(fr.colour, simple_coordinates(fr.colour), q, short):
        if not w:
            continue
        action, character, recipe = source_action(q, short, fr, mons, w)
        fixed = action - character * identity
        vectors = vectors * kernel(fixed * vectors)
        boundaries = boundaries * kernel(fixed * boundaries)
        native.append(recipe)
    return dict(
        bound=vectors.ncols() - boundaries.ncols(),
        selected_vectors=vectors.ncols(),
        invariant_boundaries=boundaries.ncols(),
        stabilizers=native,
    )


def calculate():
    t = tangent(8)
    results = []
    for degree in (6, 7):
        for record in t["orbits"]:
            q = tuple(record["root_colours"])
            for short, initial in enumerate(record["short_axis_tests"]):
                if initial["order"] is None or initial["order"] > degree:
                    continue
                first = straight(q, short, degree)
                second = curved(q, short, degree, first)
                fixed = invariant_bound(q, short, degree, first, second)
                results.append(
                    dict(
                        root_colours=q,
                        short_axis=short,
                        scalar_degree=degree,
                        straight=first,
                        curved=second,
                        finite_invariants=fixed,
                    )
                )
    return dict(weights=t, rational_tests=results)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out",
        type=Path,
        help="write deterministic gzip JSON instead of comparing to the archive",
    )
    a = p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    data = calculate()
    raw = (json.dumps(data, separators=(",", ":")) + "\n").encode()
    if a.out:
        a.out.write_bytes(gzip.compress(raw, mtime=0))
    else:
        assert raw == gzip.decompress(
            Path(__file__).with_name("nilpotent_source_bounds.json.gz").read_bytes()
        )
    for row in data["rational_tests"]:
        print(
            *row["root_colours"],
            row["short_axis"],
            row["scalar_degree"],
            row["finite_invariants"]["bound"],
        )
