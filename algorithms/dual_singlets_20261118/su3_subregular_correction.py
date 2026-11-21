"""Construct rational primitives in successive polynomial degrees of the subregular chart."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from heapq import heapify, heappop, heappush
from itertools import combinations
import json
from pathlib import Path

from su3_subregular_differential import ODD, qmonomial, variation
from su3_subregular_image import write


def quadratic_parents(word, inverse):
    parents = set()
    for i, j in combinations(range(len(word)), 2):
        key = bytes((word[i], word[j]))
        if key not in inverse:
            continue
        rest = word[:i] + word[i+1:j] + word[j+1:]
        for k in inverse[key]:
            if not ODD[k] or k not in rest:
                parents.add(bytes(sorted(rest+bytes([k]))))
    return parents


@lru_cache(6000)
def quadratic(word):
    return {v: n for v, n in qmonomial(word).items() if len(v) == len(word)+1}


def accumulate(row, column, scale):
    for key, n in column.items():
        value = row.get(key, 0) + scale*n
        if value:
            row[key] = value
        else:
            row.pop(key, None)


def reduce_row(row, pivots):
    steps = []
    if len(row) > 64:
        pending = list(row)
        heapify(pending)
        while pending:
            first = heappop(pending)
            if first not in row:
                continue
            if first not in pivots:
                break
            index, reduced = pivots[first]
            factor = row[first]
            for key, n in reduced.items():
                value = row.get(key, 0) - factor*n
                if value:
                    if key not in row:
                        heappush(pending, key)
                    row[key] = value
                else:
                    row.pop(key, None)
            steps.append((index, factor))
        return steps
    while row:
        first = min(row)
        if first not in pivots:
            break
        index, reduced = pivots[first]
        factor = row[first]
        accumulate(row, reduced, -factor)
        steps.append((index, factor))
    return steps


def solve(rhs, depth):
    inverse = defaultdict(list)
    for k in range(128):
        for pair in variation(k):
            if len(pair)==2:
                inverse[pair].append(k)
    original = dict(rhs)
    pivots, recipes = {}, []
    examined, processed = set(), set()
    frontier = set(rhs)
    for stage in range(depth):
        parents = set()
        for word in frontier:
            parents.update(quadratic_parents(word, inverse))
        parents.difference_update(processed)
        processed.update(parents)
        following = set()
        columns = []
        for parent in sorted(parents):
            row = quadratic(parent)
            following.update(row)
            if row:
                columns.append((len(row), parent, row))
        for _, parent, column in sorted(columns):
            row = dict(column)
            steps = reduce_row(row, pivots)
            if row:
                key = min(row)
                denominator = row[key]
                normalized = {t: v//denominator if isinstance(v,int) and isinstance(denominator,int)
                              and v % denominator == 0 else Fraction(v, denominator) for t,v in row.items()}
                pivots[key] = len(recipes), normalized
                recipes.append((parent, denominator, steps))
        residual = original.copy()
        steps = reduce_row(residual, pivots)
        print(stage+1, "parents", len(processed), "rank", len(pivots), "remaining", len(residual), flush=True)
        if not residual:
            coefficients = Counter(dict(steps))
            solution = {}
            for i in reversed(range(len(recipes))):
                if coefficients[i]:
                    parent, denominator, dependencies = recipes[i]
                    c = Fraction(coefficients[i], denominator)
                    solution[parent] = c
                    for j, n in dependencies:
                        coefficients[j] -= n*c
            check = original.copy()
            for parent, c in solution.items():
                accumulate(check, quadratic(parent), -c)
            assert not check
            return solution, dict(parents_examined=len(processed), independent_images=len(pivots),
                                  incidence_rounds=stage+1)
        examined.update(frontier)
        frontier = following-examined
        if not frontier:
            break
    raise ValueError("the explored quadratic images have no certified primitive")


def run(image, degree, initial, recipe, build, depth):
    raw = gzip.decompress(image.read_bytes())
    data = json.loads(raw)
    assert degree <= data["degree"]
    # Level 34 stores -Q1 O0, level 32 stores +Q1 O0. Use -Q1 O0 throughout.
    rhs = {bytes.fromhex(m):Fraction(-data["sign"]*n, data["denominator"])
           for m,n in data["coefficients"] if len(m) <= 2*degree}
    if initial:
        previous = json.loads(gzip.decompress(initial.read_bytes()))
        for encoded, c in previous:
            parent = bytes.fromhex(encoded)
            image_terms = {w:n for w,n in qmonomial(parent).items() if len(w)<=degree}
            accumulate(rhs, image_terms, -Fraction(c))
    assert all(len(w)==degree for w in rhs)
    metadata = {}
    if build:
        solution, metadata = solve(rhs, depth)
        write(recipe, [[w.hex(),str(c)] for w,c in sorted(solution.items())])
    solution = json.loads(gzip.decompress(recipe.read_bytes()))
    before = len(rhs)
    for encoded, c in solution:
        parent = bytes.fromhex(encoded)
        assert len(parent) == degree-1
        accumulate(rhs, quadratic(parent), -Fraction(c))
    assert not rhs
    return dict(level=data["level"], through_degree=degree, image_sha256=hashlib.sha256(raw).hexdigest(),
                initial_sha256=hashlib.sha256(initial.read_bytes()).hexdigest() if initial else None,
                primitive_sha256=hashlib.sha256(recipe.read_bytes()).hexdigest(), input_monomials=before,
                additional_primitive_terms=len(solution), remaining_monomials=len(rhs)), metadata


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--image", type=Path, required=True)
    p.add_argument("--degree", type=int, required=True)
    p.add_argument("--initial", type=Path)
    p.add_argument("--recipe", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--build", action="store_true")
    p.add_argument("--depth", type=int, default=3)
    args=p.parse_args()
    if not __debug__:
        p.error("verification requires assertions")
    record, metadata=run(args.image,args.degree,args.initial,args.recipe,args.build,args.depth)
    args.out.write_text(json.dumps(record,indent=2)+"\n")
    print(record, metadata, flush=True)
