"""Test six scalar insertions in the rank-one BMN threshold at both C2 root walls."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
import json
import math
from pathlib import Path
import sys

import numpy as np
from flint import nmod_mat
import mixed_incoming as m
from b2_short_image import image_samples
from dual_charge_bands import coordinates
from f4_bmn_primary import LETTERS, bmn_convention, translated
from fortuity.bmn import words as installed_words
from orthosymplectic_classical_page import minor

ALG = Path(__file__).resolve().parent.parent
PAPER = ALG / "su3_j34_20261118"
if not (PAPER/"certificate.json").is_file():
    PAPER = Path(__file__).resolve().parent
sys.path.insert(0, str(PAPER))
from verify import load_words, pairing, polynomial, polynomial_of_word
from match_threshold import paper_core

WORDS = load_words(Path(installed_words.__file__))


def insertions(even):
    values = tuple(range(0, 4, 2 if even else 1))
    for fourth in range(2):
        for second in range((6-4*fourth)//2+1):
            remaining = 6-4*fourth-2*second
            for exponents in product(values, repeat=7):
                if sum(exponents) == remaining:
                    yield fourth, second, exponents


def core_with_scalars(fourth, second, exponents):
    a,b,c,d,e,f,g = ((0,)*n for n in exponents)
    out = Counter()
    for i,j,k in permutations(range(3)):
        sign = (-1)**(int(i>j)+int(i>k)+int(j>k))
        for s,t in product(range(3), repeat=2):
            traces = ((0,0,0,0),)*fourth + ((0,0),)*second + (
                (s,)+a+(3+i,)+b, (t,)+c+(3+j,)+d,
                (3+s,)+e+(3+t,)+f+(3+k,)+g)
            factors = [WORDS.canon_word(tr) for tr in traces]
            factor_sign = math.prod(cyclic for cyclic, _ in factors)
            if factor_sign:
                phase, row = WORDS.canon_multi([tr for _, tr in factors])
                if phase:
                    out[row] += phase*factor_sign*sign
    return {row:n for row,n in out.items() if n}


@lru_cache(None)
def trace_restriction(trace, wall):
    # Long root: 2_0 + 1_1 + 1_-1. Short root: 2_1 + 2_-1.
    options = Counter()
    if wall == "long":
        sign, word = WORDS.canon_word(trace)
        if sign:
            options[0, word] += sign
        if all(letter == 0 for letter in trace) and len(trace) <= 6 and len(trace)%2 == 0:
            options[len(trace), ()] += 2
    else:
        scalar = [i for i,letter in enumerate(trace) if letter == 0]
        for n in range(0, min(6,len(scalar))+1, 2):
            for removed in combinations(scalar,n):
                remaining = tuple(x for j,x in enumerate(trace) if j not in removed)
                if len(remaining) == 1:
                    continue
                if not remaining:
                    options[n, ()] += 4
                else:
                    sign, word = WORDS.canon_word(remaining)
                    if sign:
                        options[n, word] += 2*sign
    return tuple((n,word,c) for (n,word),c in options.items() if c)


@lru_cache(maxsize=65536)
def restrict_word(word, wall):
    out = Counter()
    for choice in product(*(trace_restriction(trace,wall) for trace in word)):
        if sum(n for n,_,_ in choice) != 6:
            continue
        sign, active = WORDS.canon_multi([t for _,t,_ in choice if t])
        if sign:
            out[active] += sign*math.prod(c for _,_,c in choice)
    return tuple((w,n) for w,n in out.items() if n)


def active(expr, wall):
    result = Counter()
    for word,n in expr.items():
        for target,c in restrict_word(word,wall):
            result[target] += n*c
    return {w:n for w,n in result.items() if n}


def run(out, even, probes, workers):
    out.mkdir(exist_ok=True, parents=True)
    bmn_convention()
    certificate = json.loads((PAPER/"certificate.json").read_bytes())["channels"][2]
    ell = {(i,j):Fraction(c) for i,j,c in certificate["witness"]["monomial_coefficients"]}
    core = paper_core(WORDS)
    assert not polynomial(WORDS.Q_expr(core)) and pairing(ell,polynomial(core)) == -20
    parents = WORDS.words_of_sector((6,4,4,4),3)
    assert len(parents) == 197
    assert all(pairing(ell,polynomial(WORDS.Q_multi(w))) == 0 for w in parents)

    @lru_cache(None)
    def value(word):
        return pairing(ell,polynomial_of_word(word))

    recipes, walls, variations, seen = [], [], [], {}
    for fourth,second,powers in insertions(even):
        source = core_with_scalars(fourth,second,powers)
        row = translated(source)
        if not row:
            continue
        divisor = math.gcd(*row.values())
        if row[min(row)] < 0:
            divisor = -divisor
        key = tuple(sorted((w,n//divisor) for w,n in row.items()))
        if key in seen:
            continue
        seen[key] = len(recipes)
        recipes.append([fourth,second,powers])
        restrictions = [active(source, wall) for wall in ("long","short")]
        walls.append([sum((value(w)*n for w,n in terms.items()), Fraction())
                      for terms in restrictions])
        variations.append(m.q(row))
    print("Distinct trace polynomials",len(recipes),flush=True)
    assert any(w[0] for w in walls)
    lib = m.library(out,4,m.P)
    samples = image_samples(probes,4,820319)
    alphabet = list(LETTERS)
    matrix = np.vstack([coordinates(variations[n:n+64],4,alphabet,samples,4,"Sp",lib,workers)
                           for n in range(0,len(recipes),64)])
    residues = np.array([[w.numerator * pow(w.denominator,-1,m.P) % m.P for w in row]
                         for row in walls], dtype=np.int64)
    ranks = [
        nmod_mat(a.tolist(),m.P).rank() for a in
        (matrix, np.column_stack((matrix,residues[:,0])),
         np.column_stack((matrix,residues[:,1])), np.column_stack((matrix,residues)), residues)]
    witness = minor(matrix,m.P)
    assert witness["rank"] == ranks[0]
    wall_witness = minor(residues,m.P)
    result = dict(prime=m.P, even_insertions_only=even, recipes=recipes,
                  wall_coordinates=[[str(x) for x in row] for row in walls],
                  threshold_pairing=-20, rank_one_predecessors=len(parents),
                  image_samples=samples, tested_variation_rank=ranks[0],
                  augmented_long_rank=ranks[1], augmented_short_rank=ranks[2],
                  augmented_joint_rank=ranks[3], wall_rank=ranks[4],
                  variation_minor=witness, wall_minor=wall_witness)
    m.write(out/"test.json.gz",result)
    print("Q0, augmented long, short, both, and unconditioned wall ranks",ranks,flush=True)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--all-insertions",action="store_true")
    parser.add_argument("--probes",type=int,default=128)
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    run(args.out,not args.all_insertions,args.probes,args.workers)
