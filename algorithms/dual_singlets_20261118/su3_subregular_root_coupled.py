"""Solve the first equation modulo (f,g)^2 with all choices of its leading primitive."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path

from su3_subregular_correction import accumulate, reduce_row
from su3_subregular_differential import ODD, qmonomial, variation
from su3_subregular_image import ARCHIVES, write
from su3_subregular_root_linear import decode, root_degree


def d_variables(word):
    return bytes(k for k in word if k//16 == 5)


@lru_cache(maxsize=30000)
def projected(word):
    return {m:n for m,n in qmonomial(word).items() if root_degree(m) <= 1}


@lru_cache(None)
def inverse_rules(full=False):
    inverse = defaultdict(list)
    for k in range(128):
        for term in variation(k):
            if (full or root_degree(term) <= 1) and k//16 != 5:
                assert len(term) in (2, 3) and d_variables(term) == d_variables(bytes([k]))
                inverse[term].append(k)
    return inverse


def parents(word, image=projected, full=False):
    inverse = inverse_rules(full)
    result = set()
    for size in (2, 3):
        for positions in combinations(range(len(word)), size):
            term = bytes(word[i] for i in positions)
            if term not in inverse:
                continue
            rest = bytes(k for i,k in enumerate(word) if i not in positions)
            for k in inverse[term]:
                if ODD[k] and k in rest:
                    continue
                parent = bytes(sorted((*rest, k)))
                if (full or root_degree(parent) <= 1) and word in image(parent):
                    result.add(parent)
    return result


def serialize(row):
    return [[w.hex(),str(c)] for w,c in sorted(row.items()) if c]


def read_coefficients(rows):
    result = {bytes.fromhex(w):Fraction(c) for w,c in rows}
    if len(result) != len(rows) or not all(result.values()):
        raise ValueError("Invalid polynomial coefficient archive")
    return result


def rhs_by_spectator(path):
    data = decode(path)
    size, denominator, digest = ARCHIVES[data["level"]][1:]
    assert data["degree"] == 15 and data["linear_roots"]
    assert (data["full_image_words"], data["denominator"], data["full_image_sha256"]) == (size, denominator, digest)
    groups = defaultdict(dict)
    for w,c in data["coefficients"]:
        word = bytes.fromhex(w)
        assert root_degree(word) <= 1
        groups[d_variables(word)][word] = Fraction(-data["sign"]*c, denominator)
    return groups


def solve(rhs, image_of=projected, full=False):
    frontier = set(rhs)
    reached = set(rhs)
    examined = set()
    pivots = {}
    recipes = []
    rounds = 0
    while frontier:
        rounds += 1
        added = set()
        for word in frontier:
            added.update(parents(word, image_of, full))
        added.difference_update(examined)
        examined.update(added)
        following = set()
        for parent in sorted(added, key=lambda p:(len(image_of(p)), p)):
            image = dict(image_of(parent))
            following.update(image)
            steps = reduce_row(image, pivots)
            if image:
                first = min(image)
                scale = image[first]
                pivots[first] = len(recipes), {w:Fraction(c, scale) for w,c in image.items()}
                recipes.append((parent, scale, steps))
        residual = dict(rhs)
        steps = reduce_row(residual, pivots)
        info = dict(rounds=rounds, parent_columns=len(examined), independent_columns=len(pivots))
        if not residual:
            coefficients = Counter(dict(steps))
            solution = {}
            for i in reversed(range(len(recipes))):
                if coefficients[i]:
                    parent, scale, substitutions = recipes[i]
                    value = Fraction(coefficients[i], scale)
                    solution[parent] = value
                    for j,n in substitutions:
                        coefficients[j] -= value*n
            return dict(kind="primitive", coefficients=serialize(solution), **info)
        frontier = following-reached
        reached.update(frontier)
    key = min(residual)
    functional = {key: Fraction(1)}
    for pivot in sorted((p for p in pivots if p < key), reverse=True):
        value = -sum(n*functional.get(w, 0) for w,n in pivots[pivot][1].items() if w != pivot)
        if value:
            functional[pivot] = value
    return dict(kind="functional", coefficients=serialize(functional), **info)


def verify(rhs, result):
    vector = read_coefficients(result["coefficients"])
    if result["kind"] == "primitive":
        residual = dict(rhs)
        for word,c in vector.items():
            assert root_degree(word) <= 1
            accumulate(residual, {w:n for w,n in qmonomial(word).items() if root_degree(w)<=1}, -c)
        if residual:
            raise ValueError("The rational primitive has a nonzero residual")
        return dict(kind="primitive", terms=len(vector), residual=0)
    assert result["kind"] == "functional"
    columns = set().union(*(parents(word) for word in vector))
    if any(sum(n*vector.get(w, 0) for w,n in projected(parent).items()) for parent in columns):
        raise ValueError("The functional does not annihilate the full differential")
    pairing = sum(c*vector.get(w, 0) for w,c in rhs.items())
    if not pairing:
        raise ValueError("The functional does not separate the image")
    return dict(kind="functional", terms=len(vector), all_incident_parents=len(columns), pairing=str(pairing))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--image", type=Path, required=True)
    p.add_argument("--d-count", type=int, choices=range(5))
    p.add_argument("--spectator", help="Hexadecimal word of d jets")
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--replay", action="store_true")
    args = p.parse_args()
    if not __debug__ or (args.spectator is None) == (args.d_count is None):
        p.error("Select a spectator or d-jet count with assertions enabled")
    args.work.mkdir(parents=True, exist_ok=True)
    groups = rhs_by_spectator(args.image)
    selected = ([bytes.fromhex(args.spectator)] if args.spectator is not None else
                sorted(w for w in groups if len(w)==args.d_count))
    for word in selected:
        rhs = groups[word]
        if not rhs:
            raise ValueError("No monomials for the requested spectator")
        path = args.work/((word.hex() or "unit")+".json.gz")
        if args.replay:
            result = decode(path)
        else:
            result = solve(rhs)
        check = verify(rhs, result)
        if not args.replay:
            write(path, result)
        print(json.dumps(dict(spectator=word.hex(), input_terms=len(rhs), **check)), flush=True)
