"""Restrict the rational SO(9) source and recompute its first action at eight colours."""
import argparse
import hashlib
import math
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpq_mat, fmpz_mat

import mixed_incoming as m
from f4_pure_target import covector
from rank_four_s5 import Q1d_multi, decode_rows
from rank_four_sixth import G, IntegralFrame, pfaffian, unpack


def first(source, size):
    result = Counter()
    for word, coefficient in source.items():
        for term, factor in Q1d_multi(word, size, "SO_odd").items():
            sign, value = m.so_canon(term)
            result[value] += sign*factor*coefficient
    return {w:v for w,v in result.items() if v}


def colour_check():
    units = np.eye(64, dtype=np.int64).reshape(8,8,8,8)
    generators = np.stack([units[i,j]-units[j,i] for i in range(8) for j in range(i)])
    identity = np.eye(8, dtype=np.int64)
    lhs = -np.einsum("xij,xkl->ijkl", generators, generators)
    rhs = np.einsum("il,jk->ijkl", identity, identity)-np.einsum("ik,jl->ijkl", identity, identity)
    assert np.array_equal(lhs, rhs)
    return lhs.size


def sample(args):
    seed, words, target, sign = args
    rng = np.random.default_rng(seed)
    frame = IntegralFrame(8,1,10,[10],G,seed=1)
    for letter in G:
        matrices = np.zeros((1,10,8,8), dtype=np.int64)
        for a in matrices[0]:
            order = rng.permutation(8)
            for i,j in zip(order[::2],order[1::2]):
                direction = 2*int(rng.integers(2))-1
                a[i,j], a[j,i] = direction, -direction
        frame.L[letter] = 1, matrices
    values = frame.eval_words(words,dtype=np.int64)[:,0]
    pairing = np.einsum("iab,jba->ij",frame.L[G[0]][1][0],frame.L[G[1]][1][0])
    assert sign*int(values[target]) == math.factorial(5)*pfaffian(pairing-pairing.T)
    return values.tolist()


def run(source, primitive, workers):
    F,H = unpack(source),unpack(primitive)
    assert len(F)==11891 and len(H)==3
    sign,target = m.so_canon((G,)*5)
    boundary = m.q(H)
    assert len(boundary)==1 and target not in boundary
    scales = []
    for size in (9,8):
        image = first(F,size)
        assert set(image)==set(boundary)|{target}
        correction = next(image[w]/v for w,v in boundary.items())
        coefficient = image[target]*sign
        assert all(image[w]==correction*v for w,v in boundary.items())
        scales.append(dict(size=size,primitive_coefficient=str(correction),
                    fifth_power_coefficient=str(coefficient),
                    full_first_action=[[m.encode(w),str(v)] for w,v in sorted(image.items())]))
    assert [(r["primitive_coefficient"],r["fifth_power_coefficient"]) for r in scales]==[("1","1"),("4/5","3")]
    images = decode_rows(m.prepare((5,5,0,0,0),9),"images")
    assert all(target not in row for row in images)
    rows = images + [{target:sign}]
    words = sorted({w for row in rows for w in row})
    assert all(sum(map(len,w))==10 and len(w)<=5 for w in words)
    bound = math.factorial(10)*8**5
    assert bound<2**63
    seeds = list(range(834000,834036))
    with ProcessPoolExecutor(workers) as pool:
        values = fmpz_mat(list(pool.map(sample, [(i,words,words.index(target),sign) for i in seeds]))).transpose()
    matrix = fmpz_mat([[row.get(w,0) for w in words] for row in rows])*values
    ranks = [fmpz_mat(matrix.tolist()[:-1]).rank(),matrix.rank()]
    assert ranks[1]==ranks[0]+1
    functional = covector(matrix)
    checks = colour_check()
    return dict(group="SO(8)",reflection_parity="even",level=30,
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                primitive_sha256=hashlib.sha256(primitive.read_bytes()).hexdigest(),
                source_terms=len(F),primitive_terms=len(H),colour_completeness_entries=checks,
                actions=scales,evaluation_seeds=seeds,word_coefficient_bound=bound,
                formal_predecessors=len(images),boundary_and_augmented_ranks=ranks,
                boundary_and_target_evaluations=[[int(v) for v in row] for row in matrix.tolist()],
                separating_functional=functional)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    base=Path(__file__).with_name("rank_four")
    parser.add_argument("--source", type=Path, default=base/"source.json.gz")
    parser.add_argument("--primitive", type=Path, default=base/"primitive.json.gz")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    assert __debug__ and args.workers>0
    result=run(args.source,args.primitive,args.workers)
    m.write(args.out,result)
    print(result["boundary_and_augmented_ranks"],result["actions"][1]["fifth_power_coefficient"])
