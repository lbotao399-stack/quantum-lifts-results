"""Reconstruct the SO(9) source with fifth-power coefficient one on a fixed integer frame."""
import argparse
import math
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz, nmod_mat
from scipy.sparse import csr_matrix

import mixed_incoming as m
from mixed33_exact import read
from mixed_middle import pin
from rank_four_s5 import Q1d_multi, decode_rows
from q1d_exact import pair_terms


def rational_residue(value, modulus):
    bound = math.isqrt(modulus//2)
    a, b, c, d = modulus, 0, value, 1
    while c > bound:
        q = a//c
        a, b, c, d = c, d, a-q*c, b-q*d
    if not d or abs(d) > bound or math.gcd(c, d) != 1:
        return None
    return Fraction(c, d)


def linear(terms, weights):
    out = Counter()
    for row, c in zip(terms, weights, strict=True):
        for word, n in row.items():
            out[word] += c*n
    return {w: n for w, n in out.items() if n}


def small_matrix(rows, words, prime, workers, work):
    idx = {w: i for i, w in enumerate(words)}
    rr, cc, vv = [], [], []
    for j, row in enumerate(rows):
        for w, c in row.items():
            rr.append(j)
            cc.append(idx[w])
            vv.append(c % prime)
    coeff = csr_matrix((np.array(vv, dtype=np.int64), (rr, cc)), shape=(len(rows), len(words)))
    assert int(max(np.diff(coeff.indptr)))*(prime-1)**2 < 2**63
    alphabet = sorted({l for w in words for t in w for l in t})
    frame = m.Frame(4, 9, 10, [8], alphabet, seed=1, p=prime, group="SO_odd")
    rng = np.random.default_rng(830162)
    for letter in alphabet:
        parity, arr = frame.L[letter]
        upper = np.triu(rng.integers(-1, 2, size=arr.shape, dtype=np.int64), 1)
        frame.L[letter] = parity, (upper-upper.swapaxes(-1,-2)) % prime
    lib = m.library(work, 9, prime)
    chunks = np.array_split(np.arange(len(words)), min(2*workers, len(words)))
    with ThreadPoolExecutor(workers) as pool:
        ev = np.concatenate(list(pool.map(
            lambda js: m.values([words[j] for j in js], frame, lib), chunks)))
    test = np.linspace(0, len(words)-1, 5, dtype=int)
    np.testing.assert_array_equal(ev[test], frame.eval_words([words[j] for j in test], dtype=np.int64))
    return np.asarray(coeff @ ev) % prime


def solve_source(rows, quantum, workers, work):
    words = sorted({w for r in rows for w in r})
    assert {m.nferm(w) for w in words} == {8}
    assert {sum(map(len, w)) for w in words} == {11}
    mod, coeff, records = 1, [0]*len(rows), []
    prime = 1000003
    for _ in range(12):
        mat = small_matrix(rows, words, prime, workers, work)
        ker, dim = nmod_mat(mat.T.tolist(), prime).nullspace()
        assert dim == 1 and mat.shape == (358,405)
        v = [int(ker[i,0]) for i in range(len(rows))]
        x = sum(n*q.numerator*pow(q.denominator, -1, prime) for n,q in zip(v,quantum)) % prime
        assert x
        v = [n*pow(x,-1,prime) % prime for n in v]
        assert not (np.array(v,dtype=np.int64) @ mat % prime).any()
        record = dict(prime=prime, rank=357, matrix_sha256=m.digest(mat.astype("<u4").tobytes()))
        if not records:
            record["classical_minor"] = m.minor(mat, prime)
        records.append(record)
        inv = pow(mod,-1,prime)
        coeff = [a+mod*((b-a)*inv % prime) for a,b in zip(coeff,v)]
        mod *= prime
        candidate = [rational_residue(c,mod) for c in coeff]
        if all(c is not None for c in candidate):
            den = math.lcm(*(c.denominator for c in candidate))
            if math.gcd(den, mod) == 1 and sum(a*b for a,b in zip(candidate,quantum)) == 1:
                ints = [int(c*den) for c in candidate]
                image = linear(rows, ints)
                bound = sum(map(abs,image.values())) * math.factorial(8) * 9**11
                print("CRT", prime, "digits", len(str(mod)), "bound", bound, flush=True)
                if mod > 2*bound:
                    assert all((a*den-b) % mod == 0 for a,b in zip(coeff,ints))
                    return candidate, dict(records=records, integer_frame_seed=830162, modulus=str(mod),
                                           integer_numerator_l1=sum(map(abs,image.values())),
                                           numerator_coefficient_bound=str(bound), formal_image_words=len(image))
        prime += 2
        while not fmpz(prime).is_prime():
            prime += 2
    raise ArithmeticError("rational source did not reconstruct")


def exact_quantum(expr):
    result = Counter()
    for w,c in expr.items():
        for raw,n in Q1d_multi(w,9,"SO_odd").items():
            sign,u = m.so_canon(raw)
            if sign:
                result[u] += c*n*sign
    return {w:v for w,v in result.items() if v}


def primitive(image):
    plan = m.prepare((5,5,0,0,0),9)
    rows = decode_rows(plan, "images")
    words = sorted({w for row in rows+[image] for w in row})
    aug = fmpq_mat([[str(row.get(w,0)) for row in rows]+[str(image.get(w,0))] for w in words])
    red, rank = aug.rref()
    coeff = [Fraction(0)]*len(rows)
    for i in range(rank):
        j = next(j for j in range(red.ncols()) if red[i,j])
        assert j < len(rows)
        coeff[j] = Fraction(str(red[i,len(rows)]))
    x = linear(decode_rows(plan,"sources"),coeff)
    assert m.q(x) == image
    return x


def pack(expr):
    letters = sorted({l for w in expr for tr in w for l in tr})
    indices = {l:i+1 for i,l in enumerate(letters)}
    den = math.lcm(*(Fraction(x).denominator for x in expr.values()))
    return dict(letters=letters, denominator=den,
                words=[[bytes(i for tr in w for i in [*(indices[l] for l in tr),0]).hex(), int(v*den)]
                       for w,v in sorted(expr.items())])


def run(plan_path, first, work, workers):
    work.mkdir(parents=True, exist_ok=True)
    z, witness = read(plan_path), read(first)
    assert pin(plan_path) == witness["plan_sha256"]
    ids = witness["primes"]["1000003"]["with_first_action"]["rows"]
    rows = decode_rows(z,"images")
    q = [Fraction(witness["projected_quantum_column"][i]) for i in ids]
    coeff, checks = solve_source([rows[i] for i in ids],q,workers,work)
    recipes = decode_rows(z,"sources")
    f = linear([recipes[i] for i in ids],coeff)
    assert m.q(f) == linear([rows[i] for i in ids],coeff)
    image = exact_quantum(f)
    sign, target = m.so_canon((((1,0,()),(0,1,())),)*5)
    assert image[target] == sign
    remainder = image.copy()
    del remainder[target]
    h = primitive(remainder)
    source_letters = sorted({l for w in f for tr in w for l in tr})
    gauginos = ((1,0,()),(0,1,()))
    cross_pairs = [(l,g) for l in source_letters for g in gauginos]
    assert all(not pair_terms(l,g) and not pair_terms(g,l) for l,g in cross_pairs)
    assert all(not pair_terms(l,g) for l in gauginos for g in gauginos)
    m.write(work/"source.json.gz", pack(f))
    m.write(work/"primitive.json.gz", pack(h))
    result = dict(first_plan_sha256=pin(plan_path), source_recipe=[[i,str(c)] for i,c in zip(ids,coeff) if c],
                  replay=checks, source_terms=len(f), first_action_terms=len(image), primitive_terms=len(h),
                  source_sha256=pin(work/"source.json.gz"), primitive_sha256=pin(work/"primitive.json.gz"),
                  S5_projection="1", cross_pair_tests=2*len(cross_pairs)+len(gauginos)**2,
                  source_alphabet=source_letters)
    m.write(work/"reconstruction.json",result)
    print("source",len(f),"action",len(image),"primitive",len(h),"paired",1,flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--first", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--workers", type=int, default=8)
    a=p.parse_args()
    assert __debug__
    run(a.plan,a.first,a.work,a.workers)
