"""Test the incoming first differential on the fifth gaugino power at rank four."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from itertools import permutations
from math import ceil, comb
from pathlib import Path
import sys

import numpy as np
from flint import fmpq_mat, nmod_mat
from scipy.sparse import csr_matrix

import mixed_incoming as m
from charge_maps import content_patterns, _multiplicity, multiplicity_parts
from mixed33_exact import read
from mixed_middle import build, pin

sys.path.insert(0, str(m.RT.parent / "cluster-next-v2.nosync/work/so8-protected"))
from q1d_exact import Q1d_multi


def dimensions(group, weight, length):
    counts = Counter()
    for ds, spin_sign in ((0, 1), (1, -1)):
        for perm in permutations(range(3)):
            sign = spin_sign * (-1)**sum(i > j for n, i in enumerate(perm) for j in perm[n+1:])
            shifted = (weight[0]+ds, weight[1]-ds) + tuple(
                weight[2+i]+perm[i]-i for i in range(3))
            for bosons, fermions, count in content_patterns(shifted, length):
                mult = (_multiplicity(bosons, fermions, 0, 9) if group == "SO_odd"
                        else multiplicity_parts(bosons, fermions, 4))
                counts[sum(fermions)] += sign*count*mult
    assert all(v >= 0 for v in counts.values())
    return {d: n for d, n in sorted(counts.items()) if n}


def decode_rows(plan, kind):
    names = "source_dictionary" if kind == "sources" else "target_dictionary"
    words = list(map(m.decode, plan[names]))
    return [{words[i]: c for i, c in entries} for entries in plan[kind]]


def coordinate(w, size, target, orientation):
    image = Q1d_multi(w, size, "SO_odd")
    projected = 0
    for word, value in image.items():
        sign, canonical = m.so_canon(word)
        if sign and canonical == target:
            projected += sign*value
    return Fraction(projected, orientation)


def pure_target(prime):
    weight = (5, 5, 0, 0, 0)
    pure = decode_rows(m.prepare(weight, 10), "sources")
    lower = m.prepare(weight, 9)
    boundaries = decode_rows(lower, "images")
    sign, word = m.so_canon((((1, 0, ()), (0, 1, ())),)*5)
    assert sign and all(word not in row for row in boundaries) and not m.q({word: sign})
    words = sorted({w for row in pure+boundaries for w in row} | {word})
    def entries(rows):
        return [[r.get(w, 0) for w in words] for r in rows]
    formal_rank = fmpq_mat(entries(boundaries)).rank()
    assert len(pure) == 20 and formal_rank == 19
    result = {"formal_dimension": len(pure), "formal_boundary_rank": formal_rank,
              "orientation_in_canonical_order": sign, "groups": {}}
    for group in ("SO_odd", "Sp"):
        dims = dimensions(group, weight, 10)
        assert list(dims) == [10]
        alphabet = {l for w in words for tr in w for l in tr}
        frame = m.Frame(4, 3, 11, [10], alphabet, seed=1033405, p=prime, group=group)
        values = frame.eval_words(words, dtype=np.int64)
        base = np.asarray(entries(pure)) @ values % prime
        images = np.asarray(entries(boundaries)) @ values % prime
        target = sign*values[words.index(word)] % prime
        minor = m.minor(base, prime)
        assert minor["rank"] == dims[10]
        aug = m.minor(np.vstack((images, target)), prime)
        bd = m.minor(images, prime)
        assert bd["rank"] == formal_rank
        result["groups"][group] = dict(dimension=dims[10], cochains=minor, boundaries=bd,
                                       with_fifth_power=aug)
    return result, word, sign


def first_coordinates(plan, word, sign):
    words = list(map(m.decode, plan["source_dictionary"]))
    # A single orthogonal completeness contraction changes the number of traces
    # by at most one. Three or fewer input traces cannot produce five bilinears.
    capable = [w for w in words if len(w) >= 4]
    for n in (1, 2, 3):
        excluded = [w for w in words if len(w) == n]
        for w in excluded[:3]:
            assert coordinate(w, 9, word, sign) == 0
    coefficients = {w: coordinate(w, 9, word, sign) for w in capable}
    rows = [sum((c*coefficients.get(words[i], 0) for i, c in v), Fraction(0))
            for v in plan["sources"]]
    return rows, len(capable)


def closure_matrix(plan, prime, workers, work, seed):
    group = "SO_odd"
    weight = tuple(plan["mu"])
    dimensions_out = dimensions(group, weight, plan["length"]+1)
    assert list(dimensions_out) == [8]
    degree, width = next(iter(dimensions_out.items()))
    words = list(map(m.decode, plan["target_dictionary"]))
    assert {m.nferm(w) for w in words} == {degree}
    alphabet = {l for w in words for tr in w for l in tr}
    k = 10
    points = ceil((width+24)/comb(k, degree))
    frame = m.Frame(4, points, k, [degree], alphabet, seed=seed, p=prime, group=group)
    assert all(arr.shape[-2:] == (9, 9) for _, arr in frame.L.values())
    lib = m.library(work, 9, prime)
    chunks = np.array_split(np.arange(len(words)), min(2*workers, len(words)))
    with ThreadPoolExecutor(workers) as pool:
        values = np.concatenate(list(pool.map(
            lambda inds: m.values([words[i] for i in inds], frame, lib), chunks)))
    ids = np.linspace(0, len(words)-1, 5, dtype=int)
    np.testing.assert_array_equal(values[ids], frame.eval_words([words[i] for i in ids], dtype=np.int64))
    rr, cc, vv = [], [], []
    for j, v in enumerate(plan["images"]):
        for i, c in v:
            rr.append(j)
            cc.append(i)
            vv.append(c % prime)
    coeff = csr_matrix((np.array(vv, dtype=np.int64), (rr, cc)),
                       shape=(len(plan["images"]), len(words)))
    assert int(max(np.diff(coeff.indptr)))*(prime-1)**2 < 2**63
    return np.asarray(coeff @ values) % prime, width, dict(
        seed=seed, exterior_generators=k, points=points, evaluation_words=len(words),
        independently_evaluated_word_indices=ids.tolist())


def run(work, workers, primes, target_only):
    work.mkdir(exist_ok=True, parents=True)
    pure, word, sign = pure_target(primes[0])
    m.write(work/"pure.json", pure)
    print("pure", {k: (v["dimension"], v["boundaries"]["rank"], v["with_fifth_power"]["rank"])
                   for k, v in pure["groups"].items()}, flush=True)
    if target_only:
        return
    plan_path = build(4, work, workers, length=10, power=5)
    plan = read(plan_path)
    assert plan["mu"] == [4, 4, 1, 1, 1] and plan["source_sector"] == [1, 1, 10, 10, 10]
    col, ncapable = first_coordinates(plan, word, sign)
    result = dict(weight=plan["mu"], length=plan["length"], formal_sources=len(col),
                  source_dimensions=dimensions("SO_odd", tuple(plan["mu"]), 10),
                  plan_sha256=pin(plan_path), projected_quantum_column=list(map(str, col)),
                  evaluated_quantum_words=ncapable, nonzero_projected_rows=sum(bool(v) for v in col),
                  primes={})
    assert any(col)
    for i, prime in enumerate(primes):
        classical, dim, frame = closure_matrix(plan, prime, workers, work, 8309317+i*187)
        residual = np.array([v.numerator*pow(v.denominator, -1, prime) % prime for v in col],
                            dtype=np.int64)
        combined = np.column_stack((classical, residual))
        record = dict(classical=m.minor(classical, prime), with_first_action=m.minor(combined, prime),
                      classical_codomain_dimension=dim, frame=frame,
                      appended_column=classical.shape[1], plan_sha256=result["plan_sha256"])
        assert record["classical"]["rank"] <= dim
        assert record["with_first_action"]["rank"] <= dim+1
        result["primes"][prime] = record
        print("SO(9)", prime, dim, record["classical"]["rank"], record["with_first_action"]["rank"],
              record["with_first_action"]["determinant"], flush=True)
        m.write(work/"first.json", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--primes", nargs="+", type=int, default=[1000003, 1000033])
    p.add_argument("--target-only", action="store_true")
    args = p.parse_args()
    assert __debug__
    run(**vars(args))
