"""Cancel the ordinary scalar remainder over Q before scoring variations of a crossed source."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import gzip
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import gcd
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile
from flint import fmpz_mat

from exchange_bank import compiler, PAIRS, scalar_remainder
from exchange_scores import graph_key
from other_crossing import COLOURS, PACKET_SHA, ROOT

MIDDLES = tuple(combinations_with_replacement(range(4), 3))


def make_one(task):
    path, axes, middle = task
    sys.path.insert(0, str(path.parent))
    from build_bank import combine, dis, overlap

    raw = subprocess.run(
        ["node", str(path), json.dumps(axes), json.dumps(middle)],
        capture_output=True,
        check=True,
    )
    job = json.loads(raw.stdout)
    parts = {}
    ordinary = {}
    for side, ward in job["W"].items():
        disjoint = list(dis(job["raw"]["disjoint"]))
        parts[side] = combine(disjoint + list(overlap(ward)))
        for omitted, terms in ward.items():
            for scalar, v in scalar_remainder(terms).items():
                # Graph, orientations, ordered cycle basis and omitted ghost chords all belong to the coordinate.
                key = json.dumps(
                    [side, terms["E"], terms["Ch"], omitted, scalar],
                    separators=(",", ":"),
                )
                ordinary[key] = v
    return parts, ordinary


def solve(archive, indices, workers):
    assert sha256(archive.read_bytes()).hexdigest() == PACKET_SHA
    axes = [PAIRS[j] for j in indices]
    with TemporaryDirectory() as tmp, ZipFile(archive) as z:
        t = Path(tmp)
        source = compiler(z.read(f"{ROOT}/{COLOURS}/typed.js").decode())
        fixed = "const AX=[1,0,3,0,1,3,0,2,3,1,2,3]"
        assert source.count(fixed) == 1
        source = source.replace(
            fixed, "const AX=[1,0,3,...JSON.parse(process.argv[3]),0,2,3,1,2,3]"
        )
        (t / "generic.js").write_text(source)
        (t / "build_bank.py").write_bytes(z.read(f"{ROOT}/{COLOURS}/build_bank.py"))
        with ProcessPoolExecutor(workers) as pool:
            pieces = list(
                pool.map(make_one, [(t / "generic.js", axes, m) for m in MIDDLES])
            )
    coords = sorted(set().union(*(r for _, r in pieces)))
    if coords:
        matrix = fmpz_mat([[r.get(k, 0) for _, r in pieces] for k in coords])
        kernel, d = matrix.nullspace()
        assert d + matrix.rank() == len(MIDDLES) and (matrix * kernel).is_zero()
        vectors = [[int(kernel[i, j]) for i in range(len(MIDDLES))] for j in range(d)]
        vectors = [[x // gcd(*v) for x in v] for v in vectors]
    else:
        vectors = [
            [int(i == j) for i in range(len(MIDDLES))] for j in range(len(MIDDLES))
        ]
    bank = []
    for vec in vectors:
        sides = {}
        for side in ("left", "right"):
            weights = Counter()
            graphs = {}
            for a, (p, _) in zip(vec, pieces, strict=True):
                if a:
                    for g in p[side]:
                        k = graph_key(g)
                        weights[k] += a * g["w"]
                        graphs[k] = g
            kept = [
                {
                    **{n: graphs[k][n] for n in ("E", "d", "cp", "cm", "cubics")},
                    "w": weights[k],
                }
                for k in sorted(graphs)
                if weights[k]
            ]
            sides[side] = dict(ordinary={}, banks=kept)
        bank.append(dict(coefficients=vec, banks=sides))
    return dict(
        indices=indices,
        axes=axes,
        middle_slots=MIDDLES,
        ordinary_coordinates=len(coords),
        ordinary_rank=len(MIDDLES) - len(vectors),
        kernel=bank,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("packet", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--indices", type=int, nargs=2, required=True)
    p.add_argument("--workers", type=int, default=1)
    a = p.parse_args()
    z = solve(a.packet, a.indices, a.workers)
    a.output.write_bytes(
        gzip.compress(json.dumps(z, separators=(",", ":")).encode(), mtime=0)
    )
    print(z["indices"], z["ordinary_coordinates"], z["ordinary_rank"], len(z["kernel"]))
    print(
        [
            (len(r["banks"]["left"]["banks"]), len(r["banks"]["right"]["banks"]))
            for r in z["kernel"]
        ]
    )
