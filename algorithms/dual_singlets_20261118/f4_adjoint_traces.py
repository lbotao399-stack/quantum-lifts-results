"""Test the seven-fermion correction with traces in the 52-dimensional adjoint representation."""
import argparse
import hashlib
from functools import lru_cache
from math import lcm
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
from f4_colour import basis
from f4_first_source import multiply, singlet
from f4_formal_cycles import probe
from f4_joined_boundaries import native, parent, vary
from f4_scalar_correction import evaluate_rows
from mixed33_exact import read
from mixed_middle import pin


@lru_cache(None)
def adjoint():
    generators, _ = basis()
    gram = np.einsum("iab,jba->ij", generators, generators)
    assert np.array_equal(gram, np.diag(gram.diagonal()))
    den = lcm(*(abs(int(x)) for x in gram.diagonal()))
    scaled = []
    for x in generators:
        bracket = x@generators-generators@x
        coordinates = (np.einsum("jab,iba->ij", bracket, generators)
                       * (den//gram.diagonal())[:, None])
        assert np.array_equal(np.einsum("ij,iab->jab", coordinates, generators), den*bracket)
        scaled.append(coordinates)
    all_matrices = np.stack(scaled)
    common = int(np.gcd.reduce(np.r_[den, all_matrices.ravel()]))
    return all_matrices//common, den//common


def adjoint_frame(alphabet, points, k, degree, prime, seed):
    matrices, den = adjoint()
    generators = (matrices*pow(den, -1, prime)) % prime
    alphabet = sorted(alphabet)
    fr = m.Frame(52, points, k, [degree], alphabet, seed=1, p=prime)
    rng = np.random.default_rng(seed)
    for letter in alphabet:
        parity, shape = fr.L[letter]
        coordinates = rng.integers(0, prime, size=shape.shape[:-2]+(52,), dtype=np.int64)
        fr.L[letter] = parity, (coordinates@generators.reshape(52, -1)).reshape(shape.shape) % prime
    return fr


def run(plan, previous, work, workers):
    work.mkdir(exist_ok=True, parents=True)
    z, old = read(plan), read(previous)
    assert pin(plan) == old["formal_plan_sha256"]
    p = old["prime"]
    selected = old["ordinary_rows"]
    src_all = list(map(m.decode, z["source_dictionary"]))
    source_indices = [i for i, w in enumerate(src_all) if m.nferm(w) == 7]
    source_words = [src_all[i] for i in source_indices]
    source_lookup = {j:i for i,j in enumerate(source_indices)}
    target_words = list(map(m.decode, z["target_dictionary"]))
    target_lookup = {w:i for i,w in enumerate(target_words)}
    rhs = [(target_lookup[w], c) for w,c in m.q(multiply(singlet(), 4)).items()]
    source_rows = [[(source_lookup[j], c) for j,c in z["sources"][i]] for i in selected]
    target_rows = [z["images"][i] for i in selected]
    jordan_stage = old["stages"][1]
    libs = {n:m.library(work, n, p) for n in (27,52)}
    jordan_lib = native(work, p, "jordan")
    record = dict(prime=p, formal_plan_sha256=pin(plan),
                  previous_sha256=hashlib.sha256(previous.read_bytes()).hexdigest(),
                  adjoint_denominator=adjoint()[1], brackets_checked=52**2)
    for key, rows, words, degree, dimension, salt in (
            ("image", target_rows, target_words, 8, 144, 507100),
            ("source", source_rows, source_words, 7, 666, 507260)):
        saved = jordan_stage["ranks"][key]
        pivot_rows = saved["rows"]
        ordinary = [i for i in pivot_rows if i < len(selected)]
        seeds = [jordan_stage["seeds"][i-len(selected)] for i in pivot_rows if i >= len(selected)]
        assert pivot_rows[:len(ordinary)] == ordinary
        recipes = [rows[i] for i in ordinary] + ([rhs] if key == "image" else [])
        base, frames = probe(recipes, words, {degree:dimension}, p, libs[27], workers, salt)
        if key == "image":
            base, required = base[:-1], base[-1]
        jordan = [parent(seed, 1, "jordan", source_only=True) for seed in seeds]
        if key == "image":
            jordan = [vary(expr, 0, "jordan") for expr in jordan]
        letters = {l for w in words for tr in w for l in tr}
        extra, _ = evaluate_rows(jordan, letters, frames, p, work, workers, jordan_lib, False)
        base = np.vstack((base, extra))
        assert int(nmod_mat(base[:,saved["columns"]].tolist(), p).det()) == saved["determinant"]
        new, other_frames = probe(rows, words, {degree:dimension}, p, libs[52], workers, salt, adjoint_frame)
        assert other_frames == frames
        union = np.vstack((base, new))
        result = m.minor(union, p)
        record[key] = dict(baseline=saved["rank"], ordinary_rows=[selected[i] for i in ordinary],
                           jordan_seeds=seeds, adjoint_rows=selected, frames=frames, minor=result)
        if key == "image":
            record["with_required_image"] = m.minor(np.vstack((union, required)), p)
        m.write(work/"adjoint_traces.json", record)
        print(key, saved["rank"], result["rank"], "with RHS", record["with_required_image"]["rank"],
              flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("plan", "previous", "work"):
        parser.add_argument("--"+arg, type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    assert __debug__
    run(**vars(parser.parse_args()))
