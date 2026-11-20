"""Test the SO(7) level-24 descendant by rational colour evaluation."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from flint import fmpz_mat, nmod_mat

import mixed_incoming as m
from fortuity.full.letters import is_ferm
from gaugino_primary_order import charge
from rank_four_sixth import G
from so8_level24_descendant import raising_images, substitute
from so8_level24_generator import small_rational
from so8_instanton_sources import outside
from trace_descent import action, basis, dimensions, evaluate, highest_weights, settings


def descendant(row):
    for spin in (0, 1):
        row = substitute(row, spin, 0, False)
    singles = {tr for w in row for tr in w}
    variations = {tr: outside({(tr,): 1}, 0, True) for tr in singles}
    result = Counter()
    for w, coefficient in row.items():
        sign = 1
        for i, tr in enumerate(w):
            for (changed,), n in variations[tr].items():
                orientation, word = m.so_canon(w[:i] + (changed,) + w[i+1:])
                if orientation:
                    result[word] += orientation * sign * coefficient * n
            if sum(is_ferm(l) for l in tr) % 2:
                sign = -sign
    return {w: c for w, c in result.items() if c}


def target_frame(weight, length, seed):
    rows = highest_weights(weight, length)
    sizes = dimensions(weight, length, 7, "SO")
    letters = sorted({l for content in m.contents(weight, length) for l in content})
    _, frames = basis(rows, sizes, "SO", 7, seed, letters)
    return frames


def load_source(path):
    data = json.loads(gzip.decompress(path.read_bytes()))
    letters = [(i, j, tuple(theta)) for i, j, theta in data["letters"]]
    source = Counter()
    for numerator, denominator, traces in data["F"]:
        orientation, word = m.so_canon(tuple(tuple(letters[i] for i in tr) for tr in traces))
        source[word] += orientation * Fraction(numerator, denominator)
    return {w: c for w, c in source.items() if c}


def integral(row):
    denominator = math.lcm(*(Fraction(v).denominator for v in row.values()))
    return {w: int(v * denominator) for w, v in row.items()}, denominator


def run(source_path, seed, out):
    settings()
    source = load_source(source_path)
    parent_variation, _ = integral(m.q(source))
    parent_frame = target_frame((3,3,1,1,1), 9, seed)
    values = evaluate([parent_variation], parent_frame, 7, "SO")[0]
    assert not any(values)
    derived = descendant(source)
    closure = m.q(derived)
    assert closure == {w: -c for w, c in descendant(m.q(source)).items()}
    source_weight, length = charge(derived)
    assert (source_weight, length) == ((2,2,4,1,1), 8)
    assert all(not any(row.values()) for row in raising_images(derived))
    print("source", len(derived), "classical identity", len(closure), flush=True)

    sign, power = m.so_canon((G,)*4)
    candidate = descendant({power: sign})
    weight, length = charge(candidate)
    assert (weight, length) == ((3,3,3,0,0), 8) and not m.q(candidate)
    assert all(not any(row.values()) for row in raising_images(candidate))
    acted = action(derived, 1, 7, "SO")
    first, denominator = integral(acted)
    assert charge(acted) == (weight, length)
    print("first action", len(first), "denominator", denominator, flush=True)
    target = target_frame(weight, length, seed+6997)
    previous = highest_weights(weight, length-1)
    boundaries = [m.q(r) for rows in previous.values() for r in rows]
    coordinates = evaluate(boundaries+[candidate, first], target, 7, "SO")
    integers = [[int(v) for v in row] for row in coordinates]
    rational_ranks = [fmpz_mat(integers[:len(boundaries)+k]).rank() for k in range(3)]
    candidate_rank = rational_ranks[1]
    image_rank = fmpz_mat(integers[:-2]+integers[-1:]).rank()
    checks = []
    coefficients = []
    for prime in (1000003, 1000033):
        reduced = nmod_mat(integers[:-1], prime)
        piv = m.pivotcols(reduced)
        sub = nmod_mat([[row[j] for j in piv] for row in integers[:-1]], prime)
        indices = m.pivotcols(sub.transpose())
        assert len(piv) == candidate_rank and indices[-1] == len(boundaries)
        square = nmod_mat([[integers[i][j] for j in piv] for i in indices], prime)
        vector = nmod_mat([[integers[-1][j]] for j in piv], prime)
        solved = square.transpose().solve(vector)
        coeff = int(solved[len(indices)-1,0]) * pow(denominator, -1, prime) % prime
        coefficients.append(small_rational(coeff, prime))
        checks.append(dict(prime=prime,minor_order=len(piv),determinant=int(square.det()),
                           rows=indices,columns=piv,coefficient_mod_prime=coeff))
    assert coefficients[0] == coefficients[1] and coefficients[0]
    coefficient = coefficients[0]
    relation = [v*coefficient.denominator - denominator*coefficient.numerator*w
                for v,w in zip(integers[-1],integers[-2],strict=True)]
    relation_rank = fmpz_mat(integers[:-2]+[relation]).rank()
    assert rational_ranks == [image_rank-1, image_rank, image_rank]
    assert relation_rank == image_rank-1
    out.mkdir(parents=True, exist_ok=True)
    packed = out/"coordinates.json.gz"
    m.write(packed, dict(classical_images=integers[:-2],candidate=integers[-2],
                         first_numerator=integers[-1],first_denominator=denominator))
    record = dict(source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  colour_group="SO(7)",doubled_source_charges=m.sector(source_weight,8),
                  doubled_target_charges=m.sector(weight,8),source_grade=2,target_grade=3,
                  source_terms=len(derived),source_parent_closure_coordinates=len(values),
                  source_classical_identity_terms=len(closure),classical_closure_nonzero=0,
                  candidate_terms=len(candidate),first_action_terms=len(first),
                  first_denominator=denominator,formal_predecessors=len(boundaries),
                  target_dimension=coordinates.shape[1],rational_ranks=rational_ranks,
                  with_first_image=image_rank,with_difference=relation_rank,
                  target_coefficient=str(coefficient),finite_field_controls=checks,
                  seed=seed,frames=target,parent_frames=parent_frame,
                  coordinates_sha256=hashlib.sha256(packed.read_bytes()).hexdigest())
    m.write(out/"pairing.json", record)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=Path(__file__).resolve().parents[1]/
                        "so7_cartan_target_20261118/representatives.json.gz")
    parser.add_argument("--seed",type=int,default=633271)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    assert __debug__
    result = run(args.source, args.seed, args.out)
    print(result["rational_ranks"],result["target_coefficient"])
