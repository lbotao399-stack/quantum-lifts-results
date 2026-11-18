"""Rank all saved source profiles and solve the joint dual equation over the first prime."""
import argparse
import gzip
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from flint import nmod_mat


P = 1000003


def product(left, right):
    # Every factor is an integer in [0,p); each positive partial sum is exact.
    assert left.shape[1] == right.shape[0] and left.shape[1]*(P-1)**2 < 2**53
    return (left.astype(np.float64) @ right.astype(np.float64, copy=False)).astype(np.int64) % P


def save(path, value):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def pivots(matrix):
    reduced, rank = nmod_mat(matrix.tolist(), P).rref()
    return [next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]


class Echelon:
    def __init__(self, width, rhs=False, minors=False):
        self.width = width
        self.rhs = rhs
        self.minors = minors
        self.panels = []
        self.rank = 0
        self.witness = None

    def add(self, block, ids):
        assert block.shape == (len(ids), self.width + self.rhs)
        reduced = block.astype(np.int64)
        coefficients = []
        def subtract(panel):
            nonlocal reduced
            weight = reduced[:, panel["pivots"]].copy()
            if self.rhs:
                coefficients.append(weight)
            if weight.any():
                reduced = (reduced - product(weight, panel["rows"])) % P
        for panel in self.panels:
            subtract(panel)
        while True:
            active = np.flatnonzero(np.any(reduced[:, :self.width], axis=0))
            if not len(active):
                break
            sample = active[np.linspace(0, len(active)-1, min(len(active), len(ids)+16), dtype=int)]
            chosen = pivots(reduced[:, sample].T)
            assert chosen
            positions = sample[pivots(reduced[chosen][:, sample])]
            square = nmod_mat(reduced[chosen][:, positions].tolist(), P)
            inverse = np.asarray(square.inv().tolist(), dtype=np.int64)
            rows = product(inverse, reduced[chosen])
            assert np.array_equal(rows[:,positions], np.eye(len(chosen), dtype=np.int64))
            panel = dict(pivots=positions, rows=rows.astype(np.float64))
            if self.minors:
                panel.update(ids=[ids[i] for i in chosen], determinant=int(square.det()))
            if self.rhs:
                prior = np.column_stack([c[chosen] for c in coefficients]) if coefficients else np.zeros((len(chosen), 0), dtype=np.int64)
                panel.update(ids=[ids[i] for i in chosen], inverse=inverse.astype(np.int32), prior=prior.astype(np.int32))
            self.panels.append(panel)
            self.rank += len(positions)
            subtract(panel)
        if self.rhs and self.witness is None and reduced[:, -1].any():
            i = int(np.flatnonzero(reduced[:, -1])[0])
            scale = pow(int(reduced[i, -1]), -1, P)
            original = {ids[i]: scale}
            dual = -np.concatenate([c[i] for c in coefficients])*scale % P if coefficients else np.empty(0, dtype=np.int64)
            for panel in reversed(self.panels):
                count = len(panel["ids"])
                raw = dual[-count:] @ panel["inverse"] % P
                dual = (dual[:-count] - raw @ panel["prior"]) % P
                for tag, c in zip(panel["ids"], raw):
                    original[tag] = (original.get(tag, 0) + int(c)) % P
            assert not len(dual)
            self.witness = {k: c for k, c in original.items() if c}

    def solution(self):
        assert self.rhs and self.witness is None
        values = np.zeros(self.width, dtype=np.int64)
        assert self.width * (P - 1)**2 < 2**63
        for panel in reversed(self.panels):
            rows = panel["rows"].astype(np.int64)
            values[panel["pivots"]] = (rows[:, -1] - rows[:, :-1] @ values) % P
        return values


def controls():
    rng = np.random.default_rng(1829)
    for fail in (False, True):
        a = rng.integers(P, size=(11, 19), dtype=np.int64)
        raw = np.vstack((a[:6], a[:4], a[6:], a[:5]))
        rhs = raw @ rng.integers(P, size=19) % P
        if fail:
            rhs[-1] = (rhs[-1] + 2) % P
        joined = np.column_stack((raw, rhs))
        space = Echelon(19, True)
        for i in range(0,len(raw),5):
            space.add(joined[i:i+5], list(range(i,min(i+5,len(raw)))))
        assert space.rank == nmod_mat(raw.tolist(), P).rank()
        if fail:
            cert = np.array([space.witness.get(i,0) for i in range(len(raw))])
            assert np.array_equal(cert @ joined % P, np.r_[np.zeros(19, dtype=np.int64),1])
        else:
            assert not np.any((raw @ space.solution() - rhs) % P)


def pack(rows_path, plan, dest, workers, expected_rows, max_batches):
    dest.mkdir(parents=True, exist_ok=True)
    def pack_profile(source):
        index = source["index"]
        files = sorted(rows_path.glob("p%02d_b*.json" % index))
        assert files, index
        if expected_rows is not None:
            cap = (5*source["dimension"]+3)//4 + 32
            batches = min(max_batches, (cap+expected_rows-1)//expected_rows)
            expected = ["p%02d_b%07d" % (index, j) for j in range(batches)]
            assert [p.stem for p in files] == expected, (index, len(files), batches)
        sizes = [json.loads(path.read_text())["count"] for path in files]
        assert expected_rows is None or all(n == expected_rows for n in sizes)
        n = sum(sizes)
        headers = json.loads(files[0].read_text())["frame_widths"]
        widths = {k: val[1] for k,val in headers.items() if k != "ell_Q1"}
        arrays = {k: np.lib.format.open_memmap(dest/("p%02d_%s.npy" % (index,k)), mode="w+", dtype=np.uint32,
                                               shape=(n, w)) for k,w in widths.items()}
        rhs = np.lib.format.open_memmap(dest/("p%02d_rhs.npy" % index), mode="w+", dtype=np.uint32, shape=(n,))
        tags = []
        start = 0
        provenance = hashlib.sha256()
        for file, count in zip(files, sizes):
            raw = file.read_bytes()
            recipe = json.loads(raw)
            provenance.update(raw)
            assert recipe["index"] == index and recipe["include_source"] and file.stem=="p%02d_b%07d"%(index,recipe["serial"])
            path = file.with_suffix(".npz")
            matrix = path.read_bytes()
            assert hashlib.sha256(matrix).hexdigest() == recipe["matrix_sha256"]
            recipe_raw = file.with_suffix(".json.gz").read_bytes()
            assert hashlib.sha256(recipe_raw).hexdigest() == recipe["recipes_sha256"]
            assert len(json.loads(gzip.decompress(recipe_raw))) == count
            import io
            with np.load(io.BytesIO(matrix), allow_pickle=False) as f:
                assert set(f.files) == set(arrays) | {"ell_Q1"}
                for key, target in arrays.items():
                    value = f[key]
                    assert value.dtype == np.uint32 and value.shape == (count, widths[key]) and np.all(value < P)
                    target[start:start+count] = value
                ell = f["ell_Q1"]
                assert ell.dtype == np.uint32 and ell.shape == (count,) and np.all(ell < P)
                rhs[start:start+count] = ell
            tags.extend((recipe["serial"], j) for j in range(count))
            start += count
        assert start == n
        for arr in [*arrays.values(),rhs]:
            arr.flush()
        np.save(dest/("p%02d_tags.npy"%index),np.asarray(tags,dtype=np.int32))
        print(json.dumps(dict(packed_profile=index, rows=n)), flush=True)
        return dict(index=index, rows=n, widths=widths, receipt_sha256=provenance.hexdigest())
    with ThreadPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(pack_profile, plan["sources"]))
    save(dest/"coverage.json", records)
    return records


def run(rows, plan_path, output, batch, workers, expected_rows, max_batches):
    controls()
    plan = json.loads(plan_path.read_text())
    output.mkdir(parents=True, exist_ok=True)
    packed = output / "packed"
    records = pack(rows, plan, packed, workers, expected_rows, max_batches)
    source_ranks = {}
    closure_widths = {}
    for r in records:
        index = r["index"]
        for target, w in r["widths"].items():
            if target != "source":
                assert closure_widths.setdefault(int(target), w) == w
        source = np.load(packed/("p%02d_source.npy"%index), mmap_mode="r")
        space = Echelon(source.shape[1])
        dim = plan["sources"][index]["dimension"]
        for i in range(0, len(source), batch):
            space.add(source[i:i+batch], range(i,min(i+batch,len(source))))
            if space.rank == dim:
                break
        assert space.rank <= dim
        source_ranks[index] = space.rank
        save(output/"source_ranks.json", source_ranks)
        print(json.dumps(dict(profile=index, sampled_rank=space.rank, dimension=dim)), flush=True)
    offsets, width = {}, 0
    for t in sorted(closure_widths):
        offsets[t] = width
        width += closure_widths[t]
    space = Echelon(width, True)
    def blocks(index):
        rec = records[index]
        arrays = {int(k): np.load(packed/("p%02d_%s.npy"%(index,k)), mmap_mode="r") for k in rec["widths"] if k != "source"}
        rhs = np.load(packed/("p%02d_rhs.npy"%index), mmap_mode="r")
        for i in range(0,len(rhs),batch):
            n = min(batch,len(rhs)-i)
            mat = np.zeros((n,width+1),dtype=np.int64)
            for t, array in arrays.items():
                mat[:,offsets[t]:offsets[t]+array.shape[1]]=array[i:i+n]
            mat[:,-1]=rhs[i:i+n]
            yield i, mat
    for index in (2,0,4,5,1,9,8,7,3,6):
        for i, mat in blocks(index):
            space.add(mat, [(index,j) for j in range(i,i+len(mat))])
        save(output/"joint_progress.json", dict(latest_profile=index, closure_rank=space.rank,
                                                modular_dependency=space.witness is not None))
        print(json.dumps(dict(profile=index, joint_classical_rank=space.rank, pairing=bool(space.witness))),flush=True)
    result = dict(prime=P, job=os.environ.get("SLURM_JOB_ID"), source_rows=sum(r["rows"] for r in records),
                  source_rank=sum(source_ranks.values()), source_dimensions=[r["dimension"] for r in plan["sources"]],
                  source_ranks=source_ranks, full_source=(sum(source_ranks.values())==42991),
                  closure_rank=space.rank, augmented_rank=space.rank+int(space.witness is not None),
                  closure_frame_widths=closure_widths, eval_columns=width, field="finite field")
    if space.witness is None:
        solution=space.solution()
        np.save(output/"dual.npy",solution.astype(np.uint32))
        for index in range(10):
            for i,mat in blocks(index):
                assert not np.any((mat[:,:-1]@solution-mat[:,-1])%P)
        result["dual_residuals_checked"]=result["source_rows"]
    else:
        check=np.zeros(width+1,dtype=np.int64)
        witness=[]
        for index in range(10):
            tags=np.load(packed/("p%02d_tags.npy"%index), mmap_mode="r")
            for i,mat in blocks(index):
                factors=np.asarray([space.witness.get((index,j),0) for j in range(i,i+len(mat))], dtype=np.int64)
                check=(check+factors@mat)%P
                for local in np.flatnonzero(factors):
                    serial,sub=map(int,tags[i+local])
                    words=json.loads(gzip.decompress((rows/("p%02d_b%07d.json.gz"%(index,serial))).read_bytes()))
                    witness.append([int(factors[local]),index,serial,sub,words[sub]])
        assert np.array_equal(check,np.r_[np.zeros(width,dtype=np.int64),1])
        payload=gzip.compress(json.dumps(witness,separators=(",",":")).encode(),mtime=0)
        (output/"witness.json.gz").write_bytes(payload)
        result.update(witness_rows=len(witness), witness_sha256=hashlib.sha256(payload).hexdigest(), pairing=1)
    save(output/"result.json", result)
    print(json.dumps(result),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--rows",type=Path)
    parser.add_argument("--plan",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--batch",type=int,default=128)
    parser.add_argument("--pack-workers", type=int, default=1)
    parser.add_argument("--expected-rows-per-batch", type=int)
    parser.add_argument("--max-batches", type=int, default=5800)
    parser.add_argument("--control",action="store_true")
    args=parser.parse_args()
    if args.control:
        controls()
    else:
        run(args.rows,args.plan,args.output,args.batch,args.pack_workers,
            args.expected_rows_per_batch,args.max_batches)
