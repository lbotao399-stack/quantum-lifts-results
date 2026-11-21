"""Verify the trace factors, highest-weight recipes and Cartan determinant for omega."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import permutations, product
import json
from math import perm
from pathlib import Path

from flint import fmpz_mat

import mixed_incoming as m
from fortuity.full.letters import is_ferm
from sp3_omega_products import factor_polynomial, parity, raised_factor, raising, weight
from trace_descent import dimensions


def nonzero(terms):
    return {key: value for key, value in terms.items() if value}


def differentiate(trace, axis):
    sign = 1
    for i, (plus, minus, indices) in enumerate(trace):
        if axis < 2:
            new = (plus + (axis == 0), minus + (axis == 1), indices)
            coefficient = (plus, minus)[axis] + 1
        elif axis - 2 not in indices:
            flavour = axis - 2
            new = (plus, minus, tuple(sorted(indices + (flavour,))))
            coefficient = sign * (-1)**sum(j > flavour for j in indices)
        else:
            new, coefficient = None, 0
        if coefficient:
            yield coefficient, trace[:i] + (new,) + trace[i+1:]
        if is_ferm(trace[i]):
            sign = -sign


@lru_cache(None)
def trace_factor(spec, *, unitary=False):
    inside = tuple((0, 0, (i,)) for i, n in enumerate(spec[2:5]) for _ in range(n))
    plus, minus = (1, 0, ()), (0, 1, ())
    inside += (plus,)*spec[0] + (minus,)*spec[1]
    orders = set(permutations(inside))
    terms = {word: -1 if plus in word and minus in word and word.index(minus) < word.index(plus)
             else 1 for word in orders}
    axes = [i+2 for i in range(3) if spec[7] & (1 << i)] + [0]*spec[5] + [1]*spec[6]
    for axis in axes:
        result = Counter()
        for word, coefficient in terms.items():
            for sign, new in differentiate(word, axis):
                result[new] += coefficient * sign
        terms = nonzero(result)
    result = Counter()
    for word, coefficient in terms.items():
        if unitary:
            sign, oriented = m.fw.canon_word(word)
            tr = (oriented,)
        else:
            sign, tr = m.so_canon((word,))
        if sign:
            result[tr] += Fraction(coefficient * sign, len(orders))
    return nonzero(result)


def ordered(items, odd):
    exterior = [v for v in items if odd(v)]
    if len(set(exterior)) < len(exterior):
        return 0, ()
    inversions = sum(x > y for i, x in enumerate(exterior) for y in exterior[i+1:])
    return (-1)**inversions, tuple(sorted(items))


def odd_jet(jet):
    return not (jet[2].bit_count() % 2)


@lru_cache(None)
def diagonal(spec):
    result = Counter()
    for (trace,), coefficient in trace_factor(spec).items():
        jets = tuple((a, b, sum(1 << i for i in indices)) for a, b, indices in trace)
        sign, word = ordered(jets, odd_jet)
        if sign:
            result[word] += sign * coefficient
    result = nonzero(result)
    assert result == factor_polynomial(spec)
    assert all(n.denominator == 1 for n in result.values())
    return {word: int(n) for word, n in result.items()}


def check_factor(spec):
    factor = trace_factor(spec)
    assert not m.q(factor)
    for axis, simple in enumerate((3, 1, 2)):
        difference = Counter()
        for tr, coefficient in factor.items():
            for target, n in m.raise_word(tr, simple).items():
                difference[target] += coefficient*n
        for n, raised in raised_factor(spec, axis):
            for target, coefficient in trace_factor(raised).items():
                difference[target] -= n*coefficient
        assert not nonzero(difference), (spec, axis)
    diagonal(spec)


@lru_cache(None)
def partitions(length):
    # Set partitions in order of first occurrence of each colour.
    choices = []
    for labels in product(range(3), repeat=length):
        if list(dict.fromkeys(labels)) == list(range(max(labels)+1)):
            choices.append((labels, tuple(tuple(i for i, v in enumerate(labels) if v == block)
                                           for block in range(max(labels)+1))))
    return choices


@lru_cache(maxsize=12000)
def one_colour(factors):
    if not factors:
        return {(): 1}
    if len(factors) == 1:
        return diagonal(factors[0])
    values = Counter()
    for w, n in one_colour(factors[:-1]).items():
        for v, c in diagonal(factors[-1]).items():
            sign, combined = ordered(w + v, odd_jet)
            if sign:
                values[combined] += n*c*sign
    return nonzero(values)


@lru_cache(None)
def odd_slot(slot):
    return sum(odd_jet(jet) for jet in slot) % 2


def product_coordinates(factors, columns):
    result = Counter()
    for labels, blocks in partitions(len(factors)):
        assignments = perm(3, len(blocks))
        odd = [labels[i] for i, factor in enumerate(factors) if parity(factor)]
        swaps = sum(x > y for i, x in enumerate(odd) for y in odd[i+1:])
        scale = (-1)**swaps * assignments * 2**len(factors)
        slots = [one_colour(tuple(factors[i] for i in block)).items() for block in blocks]
        for monomials in product(*slots):
            words = tuple(w for w, _ in monomials) + ((),)*(3-len(blocks))
            sign, column = ordered(words, odd_slot)
            if sign and column in columns:
                coefficient = scale*sign
                for _, n in monomials:
                    coefficient *= n
                result[columns[column]] += coefficient
    return nonzero(result)


def run(recipe_path, output):
    raw = recipe_path.read_bytes()
    data = json.loads(gzip.decompress(raw))
    assert data["group"] == "USp(6)"
    recipes = [[(tuple(map(tuple, w)), n) for w, n in row] for row in data["highest_weight_recipes"]]
    columns = {tuple(tuple(map(tuple, s)) for s in col): i for i, col in enumerate(data["cartan_columns"])}
    assert len(recipes) == len(columns) == data["cartan_rank"] > 0
    products = sorted({w for row in recipes for w, _ in row})
    factors = sorted({g for w in products for g in w})
    for spec in factors:
        check_factor(spec)
    print("Closed factors and raising identities", len(factors), flush=True)
    coordinates = {}
    for w in products:
        assert sum(sum(g[:5]) for g in w) == data["length"]
        assert [sum(values) for values in zip(*(weight(g) for g in w))] == data["weight"]
        coordinates[w] = product_coordinates(w, columns)
    matrix = []
    for row in recipes:
        highest, values = Counter(), Counter()
        for w, n in row:
            for key, coefficient in raising(w).items():
                highest[key] += n*coefficient
            for j, coefficient in coordinates[w].items():
                values[j] += n*coefficient
        assert not nonzero(highest)
        matrix.append([values[j] for j in range(len(columns))])
    determinant = int(fmpz_mat(matrix).det())
    witness = data["minor"]
    assert witness["columns"] == witness["rows"] == list(range(len(columns)))
    assert determinant % data["prime"] == witness["determinant"] != 0
    source_character = dimensions(tuple(data["weight"]), data["length"], 6, "Sp")
    result = dict(recipe_sha256=hashlib.sha256(raw).hexdigest(), closed_trace_factors=len(factors),
                  verified_factor_raising_equations=3*len(factors), products=len(products),
                  highest_weight_recipes=len(recipes), source_character=source_character,
                  cartan_rank=len(columns), determinant=str(determinant),
                  determinant_residues={str(p): determinant % p for p in (1000003, 1000033)})
    m.write(output, result)
    print("Cartan rank and residues", len(columns), result["determinant_residues"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert __debug__
    run(args.products, args.out)
