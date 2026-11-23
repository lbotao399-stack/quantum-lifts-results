"""Check the two exterior signs, the finite-field sections, and the actual determinants."""

import argparse
import ast
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import tempfile
import zipfile

import numpy as np
from flint import nmod_mat

ARCHIVE_SHA256 = "b668a97232466286970c3c3b7df49128d5b98e8e0a4b00e23f57fa75563f7972"
DIRECTORY = "f4-j36-eight-Albert-certificate"
GRAPHS = "f4-j36-eight-cubic"
PIVOTS = "f4-j36-nested1-seal"
POLE_WEIGHTS = (-4 * 2392 * 7113106, 2 * 2299 * 1850212)
DENOMINATOR = 2 * Fraction(28, 3) ** 4
TRACE_SCORES = (
    -104166480711910256640000,
    -17361080118651709440000,
    -109953507418127493120000,
    0,
)
NESTED_SCORE = Fraction(-598031033628600, 7)
POLE_LINE = re.compile(r"graph (\d+) poles \[(\d+), (\d+)\] weighted raw (\d+)")


def read(path):
    return json.loads(path.read_bytes())


def mod(x, prime):
    return x.numerator * pow(x.denominator, -1, prime) % prime


def leaves(tree):
    if isinstance(tree, int):
        yield tree
    else:
        for child in tree:
            yield from leaves(child)


def parity(order):
    return (-1) ** sum(a > b for i, a in enumerate(order) for b in order[i + 1 :])


def exterior_sign(graph, owners, tree):
    edges = sorted(set(graph["cp"]) | set(graph["cm"]))
    owner = {e: graph["E"][e][side] for e, side in zip(edges, owners, strict=True)}
    traversal = [e for vertex in leaves(tree) for e in edges if owner[e] == vertex]
    assert sorted(traversal) == edges
    sign = 1
    for declared in (graph["cp"], graph["cm"]):
        sign *= parity([declared.index(e) for e in traversal if e in declared])
    return sign


def records(path):
    plans = {}
    for line in path.read_text().splitlines():
        if line.startswith("plan "):
            _, index, tail = line.split(" ", 2)
            _, owners, tree = ast.literal_eval(tail.rsplit(" ", 1)[0])
            plans[int(index)] = owners, tree
        match = POLE_LINE.match(line)
        if match:
            index, a, b, weighted = map(int, match.groups())
            owners, tree = plans[index]
            yield index, a, b, weighted, owners, tree


def graph_totals(folder, banks, side, prime):
    graphs = banks[side]["banks"]
    source_sides = ["left", "right"] if side == "right" else ["left"]
    values = {}
    for src in source_sides:
        for path in sorted(folder.glob(f"native-B-{src}-*-{prime}.log")):
            for i, a, b, reported, owners, tree in records(path):
                if src != side and i >= 12:
                    continue
                g = graphs[i]
                if src != side:
                    assert all(
                        g[k] == banks[src]["banks"][i][k]
                        for k in ("E", "d", "w", "cp", "cm")
                    )
                raw = (POLE_WEIGHTS[0] * a + POLE_WEIGHTS[1] * b) * g["w"] % prime
                assert raw == reported
                pair = (raw, raw * exterior_sign(g, owners, tree) % prime)
                assert i not in values or values[i] == pair
                values[i] = pair
    # This reproduces the upstream mirror rule. The optional native run below
    # evaluates every mirror directly at characteristic zero.
    by_pair = {
        (tuple(g["cp"]), tuple(g["cm"]), str(g["E"]), g["w"]): i
        for i, g in enumerate(graphs)
    }
    for i, g in enumerate(graphs):
        mirror = by_pair[tuple(g["cm"]), tuple(g["cp"]), str(g["E"]), g["w"]]
        if mirror in values:
            assert i not in values or values[i] == values[mirror]
            values[i] = values[mirror]
    assert len(values) == len(graphs)
    return [sum(pair[k] for pair in values.values()) % prime for k in (0, 1)]


def section_and_matrix(root, prime):
    pivot = root / PIVOTS
    graph = root / GRAPHS
    chart = read(pivot / "selected-1009.json")
    rows = chart["rows"]
    nested = rows.index(["nested", 0])
    raw = gzip.decompress((pivot / f"selected-{prime}.u16.gz").read_bytes())
    a = nmod_mat(
        np.frombuffer(raw, dtype="<u2").reshape(703, 703).astype(int).tolist(), prime
    )
    y = [row[1] for row in read(graph / f"outgoing-{prime}.json")]
    stored = read(graph / f"coordinates-{prime}.json")["coordinates"][1]
    coeff = nmod_mat([stored], prime)
    assert coeff * a == nmod_mat([y], prime)
    scores = [
        mod(
            sum(
                Fraction(value) * s for value, s in zip(row, TRACE_SCORES, strict=True)
            ),
            prime,
        )
        for row in read(pivot / "trace_C4_ordered.json")["rows"]
    ]
    scores[nested] = mod(NESTED_SCORE, prime)
    section = sum(v * s for v, s in zip(stored, scores, strict=True)) % prime
    return a, y, scores, section


def augmented(a, y, s, value, p):
    rows = a.tolist()
    for row, score in zip(rows, s, strict=True):
        row.append(score)
    return int(nmod_mat(rows + [y + [value]], p).det())


def check(root):
    banks = read(root / GRAPHS / "bank-B-reduced.json")
    result = []
    for p in (1009, 1013):
        a, y, s, section = section_and_matrix(root, p)
        for side in ("left", "right") if p == 1009 else ("left",):
            raw, ordered = graph_totals(root / GRAPHS, banks, side, p)
            old, correct = (v * mod(1 / DENOMINATOR, p) % p for v in (raw, ordered))
            before, after = (augmented(a, y, s, v, p) for v in (old, correct))
            assert int(a.det()) == {1009: 1005, 1013: 734}[p]
            assert (
                before == {1009: 551, 1013: 478}[p]
                and correct == section
                and after == 0
            )
            result.append(
                dict(
                    prime=p,
                    mark=side,
                    classical_det=int(a.det()),
                    archived_raw=raw,
                    corrected_raw=ordered,
                    archived_score=old,
                    corrected_score=correct,
                    section=section,
                    archived_augmented_det=before,
                    corrected_augmented_det=after,
                )
            )
    return result


def native_replay(root):
    # Use the pinned original evaluator, with an independent order conversion.
    import os
    import sys

    os.environ["P"] = "0"
    folder = root / GRAPHS
    sys.path.insert(0, str(folder))
    import native

    graph_banks = read(folder / "bank-B-reduced.json")
    plans = {}
    for side in ("left", "right"):
        for path in sorted(folder.glob(f"native-B-{side}-*-1009.log")):
            for i, _, _, _, bits, tree in records(path):
                g = graph_banks[side]["banks"][i]
                plans.setdefault(
                    (str(g["E"]), tuple(sorted((tuple(g["cp"]), tuple(g["cm"]))))),
                    (bits, tree),
                )

    def evaluate(g, owners, tree):
        previous = native.plan
        try:
            native.plan = lambda *unused: (0, owners, tree)
            poles = [native.run(g, j) for j in (0, 1)]
        finally:
            native.plan = previous
        sign = exterior_sign(g, owners, tree)
        return poles, sign * g["w"] * sum(
            x * y for x, y in zip(POLE_WEIGHTS, poles, strict=True)
        )

    graph_scores = {}
    totals = {}
    for side in ("left", "right"):
        all_scores = []
        for i, g in enumerate(graph_banks[side]["banks"]):
            key = str((g["E"], g["d"], g["cp"], g["cm"], g["w"]))
            if key not in graph_scores:
                bits, tree = plans[
                    str(g["E"]), tuple(sorted((tuple(g["cp"]), tuple(g["cm"]))))
                ]
                poles, score = evaluate(g, bits, tree)
                if (side, i) == ("left", 0):
                    order = (tree[0], (tree[1][1], tree[1][0]))
                    other, check_score = evaluate(g, bits, order)
                    # Any plan with this exact tree is sufficient for the swap control.
                    assert (
                        poles == [-5750, 24200]
                        and other == [5750, -24200]
                        and score == check_score
                    )
                graph_scores[key] = score
            all_scores.append(graph_scores[key])
        totals[side] = sum(all_scores) / DENOMINATOR
    assert totals["left"] == totals["right"] == Fraction(-299015516814300, 49)
    return {k: str(v) for k, v in totals.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("zip", type=Path)
    parser.add_argument(
        "--native",
        action="store_true",
        help="regenerate both complete graph banks over the rationals",
    )
    args = parser.parse_args()
    if not __debug__ or sha256(args.zip.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise SystemExit("use the pinned archive and enable Python assertions")
    with tempfile.TemporaryDirectory() as temp, zipfile.ZipFile(args.zip) as z:
        z.extractall(temp)
        root = Path(temp) / DIRECTORY
        result = {"modular": check(root)}
        if args.native:
            result["exact_graph_scores"] = native_replay(root)
        print(json.dumps(result, indent=2))
