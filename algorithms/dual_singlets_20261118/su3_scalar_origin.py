"""Reduce the leading scalar degree of the signed SU(3) first image on unrestricted matrices."""
import argparse
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import sys

from su3_subregular_image import ALGORITHMS, ARCHIVES, write

sys.path.insert(0, str(ALGORITHMS / "su3_j34_20261118"))
import formal_trace_blocks as formal


SCALARS = {formal.ID[0, 0, (i,)] for i in range(3)}
GAUGINOS = {formal.ID[0, 1, ()], formal.ID[1, 0, ()]}


@lru_cache(maxsize=100000)
def scalar_count(trace):
    return sum(k in SCALARS for k in trace)


def central_trace(trace):
    if not any(k in GAUGINOS for k in trace):
        return 1, trace
    ghosts = bytes(sorted(k for k in trace if k in GAUGINOS))
    body = bytes(k for k in trace if k not in GAUGINOS)
    inversions = sum(formal.ODD[k] for i, g in enumerate(trace) if g in GAUGINOS for k in trace[:i]
                     if k not in GAUGINOS or k > g)
    sign = (-1)**inversions
    if len(body) > 1:
        phase, body = formal.canon(body)
        if not phase:
            return 0, b""
        sign *= phase
    phase, trace = formal.canon(ghosts + body)
    return sign * phase, trace


def add(values, word, coefficient):
    value = values.get(word, 0) + coefficient
    if value:
        values[word] = value
    else:
        values.pop(word, None)


@lru_cache(maxsize=50000)
def relations(trace):
    count = scalar_count(trace)
    parents = set()
    for i in range(len(trace)):
        rotated = trace[i:] + trace[:i]
        if not any(k in GAUGINOS for k in rotated[:2]):
            continue
        merged = formal.MERGE.get(rotated[:2])
        if merged is None:
            continue
        sign, parent = formal.canon(bytes((merged,)) + rotated[2:])
        if sign:
            parents.add(parent)
    out = {}
    for parent in sorted(parents):
        image = {t: c for c, t in formal.qtrace(parent)}
        kept = {t: n for t, n in image.items() if scalar_count(t) <= count and n}
        if trace not in kept:
            continue
        if len(kept) > 2 or any(scalar_count(t) != count for t in kept):
            raise ArithmeticError("The leading variation is not a single graded interchange")
        out[parent] = kept
    return out


def straightener():
    reductions = {}

    def reduce(trace):
        if trace in reductions:
            return reductions[trace]
        if not any(k in GAUGINOS for k in trace):
            return 1, trace, {}
        expected, canonical = central_trace(trace)
        root = canonical if expected else trace
        # Each trace is its multiple of root plus the differential of its path in the interchange graph.
        paths = {root: (1, {})}
        root_primitive = None
        queue = deque([root])
        while queue:
            current = queue.popleft()
            weight, primitive = paths[current]
            for parent, row in relations(current).items():
                new = next((t for t in row if t != current and t not in paths), None)
                if new is not None:
                    factor = Fraction(-row[current], row[new])
                    recipe = {k: factor * v for k, v in primitive.items()}
                    add(recipe, parent, Fraction(1, row[new]))
                    paths[new] = factor*weight, recipe
                    queue.append(new)
                    continue
                if not expected and root_primitive is None:
                    scale = sum(c * paths[t][0] for t, c in row.items())
                    if scale:
                        recipe = {parent: Fraction(1, scale)}
                        for t, c in row.items():
                            for k, v in paths[t][1].items():
                                add(recipe, k, -Fraction(c, scale)*v)
                        root_primitive = recipe
        if trace not in paths or (not expected and root_primitive is None):
            raise ArithmeticError("The cyclic normal form is not connected by the recorded variations")
        for t, (weight, recipe) in paths.items():
            if root_primitive is not None:
                for k, v in root_primitive.items():
                    add(recipe, k, weight*v)
                weight = 0
            sign, normal = central_trace(t)
            if weight != sign or (weight and normal != root):
                raise ArithmeticError("The differential and the graded permutation give different signs")
            reductions[t] = weight, normal, recipe
        return reductions[trace]
    return reduce


@lru_cache(maxsize=100000)
def encode_trace(trace):
    return bytes(formal.ID[letter] for letter in trace)


def leading(level, all_terms=False):
    relative, size, denominator, digest = ARCHIVES[level]
    degree = {32: 5, 34: 4}[level]
    counts = Counter()
    terms = {}
    sha = hashlib.sha256()
    with gzip.open(ALGORITHMS / relative, "rb") as stream:
        for line in stream:
            sha.update(line)
            coefficient, w = json.loads(line)
            word = tuple(encode_trace(tuple((z, t, tuple(theta)) for z, t, theta in tr)) for tr in w)
            n = sum(map(scalar_count, word))
            counts[n] += 1
            if n != degree and not all_terms:
                continue
            if n == degree and sum(k in GAUGINOS for tr in word for k in tr) != 2:
                raise ValueError("Unexpected derivative content at the lowest scalar degree")
            value = (1 if level == 32 else -1) * Fraction(coefficient) * denominator
            if word in terms or value.denominator != 1 or not value:
                raise ValueError("Invalid trace coefficient")
            terms[word] = value.numerator
    if sha.hexdigest() != digest or counts.total() != size or min(counts) != degree:
        raise ValueError("Source digest, term count or minimum scalar degree disagrees")
    return degree, counts, denominator, terms


def reduce_image(terms):
    reduce = straightener()
    correction = {}
    remaining = {}
    for word, value in terms.items():
        factors = list(word)
        prefix = 1
        for i, trace in enumerate(word):
            weight, normal, recipe = reduce(trace)
            for parent, c in recipe.items():
                factors[i] = parent
                sign, joined = formal.multi(tuple(factors))
                if sign:
                    add(correction, joined, -value*prefix*sign*c)
            if not weight:
                break
            factors[i] = normal
            prefix *= (-1)**formal.par(normal)
            value *= weight
        else:
            sign, joined = formal.multi(tuple(factors))
            if sign:
                add(remaining, joined, sign*value)
    return correction, remaining


def serialize(poly, denominator):
    return [[[t.hex() for t in w], str(Fraction(n, denominator))] for w, n in sorted(poly.items())]


def verify(terms, correction, degree):
    residual = dict(terms)
    for word, value in correction.items():
        for target, c in formal.q(word).items():
            if sum(map(scalar_count, target)) <= degree:
                add(residual, target, c*value)
    return residual


def read(poly, denominator):
    result = {}
    for strands, coefficient in poly:
        word = tuple(bytes.fromhex(t) for t in strands)
        value = Fraction(coefficient)*denominator
        if not value or word in result:
            raise ValueError("Invalid correction archive")
        result[word] = value
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, choices=(32, 34), required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    args.data.mkdir(exist_ok=True, parents=True)
    degree, counts, den, terms = leading(args.level)
    archive = args.data / f"j{args.level}_leading.json.gz"
    if args.replay:
        saved = json.loads(gzip.decompress(archive.read_bytes()))
        correction, normal = (read(saved[name], den) for name in ("correction", "residual"))
    else:
        correction, normal = reduce_image(terms)
    residual = verify(terms, correction, degree)
    if residual != normal:
        raise ArithmeticError("The original graded trace differential disagrees with the straightening")
    if not args.replay:
        write(archive, dict(correction=serialize(correction, den), residual=serialize(residual, den)))
    info = dict(level=args.level, source_sha256=ARCHIVES[args.level][3], input_by_scalar_degree=dict(sorted(counts.items())),
                lowest_scalar_degree=degree, leading_input=len(terms), correction_terms=len(correction),
                residual_terms=len(residual), rational_replay=True)
    path = args.data / f"j{args.level}_leading.json"
    if args.replay:
        if json.loads(path.read_bytes()) != json.loads(json.dumps(info)):
            raise ValueError("Receipt disagrees")
    else:
        path.write_text(json.dumps(info, indent=2) + "\n")
    print(info)
