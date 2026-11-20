"""Rank the incoming classical differential in the two remaining mixed singlet sectors.

Each trace shape is an invariant subspace of all three spin/flavour raising
operators. Integer dependencies give a rational kernel basis. The output
contains the original rational source recipes and their differential, so
each lower rank bound can be replayed directly in finite colour.
"""
import argparse
import gzip
import hashlib
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mixed_incoming as m
from mixed33_exact import apply, exact


def pin(path):
    with gzip.open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def solve_shape(task):
    index, words, path = task
    if path.exists():
        return index
    equations = [{(axis, term): c for axis in (1, 2, 3)
                  for term, c in m.raise_word(w, axis).items()} for w in words]
    keys = {v: i for i, v in enumerate(sorted({v for eq in equations for v in eq}))}
    rows = [{keys[v]: c for v, c in eq.items()} for eq in equations]
    independent, _, kernel = exact(rows)
    assert len(independent) + len(kernel) == len(words)
    sources, images = [], []
    for vector in kernel:
        assert not apply(vector, rows)
        source = {words[i]: c for i, c in vector.items()}
        image = m.q(source)
        assert not m.q(image)
        sources.append(source)
        images.append(image)
    sw, sr = m.flatten(sources)
    tw, tr = m.flatten(images)
    m.write(path, dict(index=index, words=len(words), shape=m.shape(words[0], "spin"),
                       dimension=len(kernel), source_dictionary=list(map(m.encode, sw)),
                       sources=sr, target_dictionary=list(map(m.encode, tw)), images=tr))
    return index


def build(t, work, workers, length=None):
    if length is None:
        length = 10 - t
    mu = (t, t, 4-t, 4-t, 4-t)
    assert m.sector(mu, 11-t) == (1, 1, 8, 8, 8)
    plan = work / "plan.json.gz"
    if plan.exists():
        return plan
    blocks = m.words_by_shape(length, mu)
    ordered = [v for k, v in sorted(blocks.items())]
    total = sum(map(len, ordered))
    print("source", mu, length, "words", total, "shapes", len(ordered),
          "largest", max(map(len, ordered)), flush=True)
    paths = [work / f"shape-{i:05}.json.gz" for i in range(len(ordered))]
    tasks = sorted(((i, v, paths[i]) for i, v in enumerate(ordered)),
                   key=lambda r: -len(r[1]))
    with ProcessPoolExecutor(workers) as pool:
        for _ in pool.map(solve_shape, tasks, chunksize=1):
            pass
    sources, images = [], []
    for i, (original, path) in enumerate(zip(ordered, paths)):
        z = json.loads(gzip.decompress(path.read_bytes()))
        assert z["index"] == i and z["words"] == len(original)
        assert z["dimension"] == len(z["sources"]) == len(z["images"])
        for key, vec, out in (("source_dictionary", "sources", sources),
                              ("target_dictionary", "images", images)):
            ws = list(map(m.decode, z[key]))
            out.extend({ws[k]: c for k, c in row} for row in z[vec])
    sw, sr = m.flatten(sources)
    tw, tr = m.flatten(images)
    m.write(plan, dict(mu=mu, length=length, source_words=total, source_shapes=len(paths),
                       formal_singlets=len(sources), source_sector=m.sector(mu, length),
                       target_sector=m.sector(mu, length+1), source_dictionary=list(map(m.encode, sw)),
                       sources=sr, target_dictionary=list(map(m.encode, tw)), images=tr))
    print("formal singlets", len(sources), "image words", len(tw), flush=True)
    return plan


def run(t, work, workers, primes, dimensions):
    work.mkdir(exist_ok=True, parents=True)
    plan = build(t, work, workers)
    doc = json.loads(gzip.decompress(plan.read_bytes()))
    digest = pin(plan)
    by_group = {r["group"]: r for r in json.loads(dimensions.read_bytes()) if r["mu"] == doc["mu"]}
    for prime in primes:
        for group in ("SO_odd", "Sp"):
            out = work / f"{group}-{prime}.json"
            if out.exists():
                assert json.loads(out.read_bytes())["plan_sha256"] == digest
                continue
            target = by_group[group]["by_length_and_fermion_degree"][str(11-t)]
            degrees = sorted(set(m.nferm(m.decode(w)) for w in doc["target_dictionary"]))
            assert set(degrees) == set(map(int, target))
            k = 8
            points = [math.ceil((target[str(d)] + 24)/math.comb(k, d)) for d in degrees]
            frame_dir = work / f"frame-{group}-{prime}"
            frame_dir.mkdir(exist_ok=True)
            z = m.run(doc, group, prime, frame_dir, workers, points, k,
                      7311322 + t*97 + prime, 3)
            z.update(plan_sha256=digest, mu=doc["mu"], source_length=doc["length"],
                     finite_colour_dimension=sum(target.values()))
            m.write(out, z)
            print("finite-colour rank", t, group, prime, z["rank"], z["determinant"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", choices=[1, 2, 3], type=int, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--primes", type=int, nargs="*", default=[1000003, 1000033])
    parser.add_argument("--dimensions", type=Path, default=Path(__file__).with_name("charge_maps.json"))
    assert __debug__
    run(**vars(parser.parse_args()))
