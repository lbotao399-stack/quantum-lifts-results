"""Rebuild the native finite-colour matrices and their selected determinants."""
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
from flint import nmod_mat

from build_g2_full_chain import build
from g2_633_relation import convert


HERE = Path(__file__).resolve().parent


def pivots(matrix, prime):
    echelon, rank = nmod_mat(matrix.tolist(), prime).rref()
    columns = []
    column = 0
    for row in range(rank):
        while not echelon[row, column]:
            column += 1
        columns.append(column)
        column += 1
    return columns


def certify(engine, packet, weight, y, prime, seed, workers, scratch):
    raw = scratch / "matrix.u32"
    env = dict(os.environ, G2_PRODUCT_MATRIX=str(raw))
    cmd = [engine, *weight, y, seed, workers, prime, packet]
    record = json.loads(subprocess.check_output(list(map(str, cmd)), env=env))
    data = np.fromfile(raw, dtype="=u4")
    p, nrows, nbound, nprod = map(int, data[:4])
    assert p == prime and record["weight"] == list(weight[3:]+weight[:3])
    assert (record["dimension"], record["products"]) == (nrows, nprod)
    matrix = data[4:].reshape(nrows, nbound+nprod)
    col = pivots(matrix, prime)
    assert [i-nbound for i in col if i >= nbound] == record["selected"]
    for field, subset in (("boundary_minor", [i for i in col if i < nbound]), ("joined_minor", col)):
        row = pivots(matrix[:, subset].T, prime)
        assert len(row) == len(subset)
        square = np.ascontiguousarray(matrix[np.ix_(row, subset)], dtype="<u4")
        det = int(nmod_mat(square.tolist(), prime).det())
        assert det and len(row) == record["boundary" if field == "boundary_minor" else "augmented"]
        record[field] = dict(rows=row, columns=subset, det=det, sha256=sha256(square.tobytes()).hexdigest())
    record["predecessor"] = nbound
    return record


def regenerate(repo, workers):
    def run(*args):
        subprocess.run(list(map(str, args)), check=True)

    with TemporaryDirectory() as tmp:
        dest = Path(tmp)
        build(repo, dest)
        engine = dest / "product"
        run("c++", "-std=c++17", "-O3", "-pthread", "-I" + str(dest), HERE/"g2_product_span.cpp", "-o", engine)
        alg = repo / "algorithms/dual_singlets_20261118"
        convert(alg / "taylor_hooks/bmn_633", dest/"P0.dat")
        tasks = [
            ("basis", (6,3,3,0,0), 10, "g2-633-product-recipes.json.gz"),
            ("partner-t1", (5,2,2,1,1), 8, "g2-522-recipes.json.gz"),
            ("partner-t3", (3,0,0,3,3), 6, "g2-633-partner-t3-recipes.json.gz"),
        ]
        for key, weight, y, recipes in tasks:
            run(sys.executable, HERE/"g2_product_words.py", "--repo", repo, "--recipes", HERE/recipes,
                "--out", dest/(key+".dat"))
        plans = [
            ("basis-2147483647", "basis", (6,3,3,0,0), 10, 2147483647,883149),
            ("basis-1000033", "basis", (6,3,3,0,0), 10, 1000033,701773),
            ("P0-2147483647", "P0", (6,3,3,0,0), 10, 2147483647,881753),
            ("P0-1000033", "P0", (6,3,3,0,0), 10, 1000033,700169),
            ("partner-t1", "partner-t1", (5,2,2,1,1), 8, 2147483647,691447),
            ("partner-t3", "partner-t3", (3,0,0,3,3), 6, 2147483647,953279),
        ]
        output = {}
        for name, key, w, y, p, seed in plans:
            print("Checking", name, flush=True, file=sys.stderr)
            output[name] = certify(engine, dest/(key+".dat"), w, y, p, seed, workers, dest)
        ranks = json.loads((HERE/"partner-ranks.json").read_bytes())
        for t, seed in enumerate((731099, 934189, 953089, 953251)):
            w = (6-t,3-t,3-t,t,t)
            length = 9-t
            cmd = [dest/"g2_dotted_chain",*w[:3],length,length,seed,workers,2147483647,"--probe-last",*w[3:]]
            result = json.loads(subprocess.check_output(list(map(str,cmd))))
            block, = result["blocks"]
            assert result["weight"] == list(w) and block["length"] == length
            assert [block[key] for key in ("dimension","incoming","outgoing")] == ranks[str(t)]
        return output


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=2)
    a = p.parse_args()
    if not __debug__ or a.workers < 1:
        raise ValueError("enable assertions and select positive worker count")
    result = regenerate(a.repo.resolve(), a.workers)
    a.out.write_text(json.dumps(result,indent=2)+"\n")
