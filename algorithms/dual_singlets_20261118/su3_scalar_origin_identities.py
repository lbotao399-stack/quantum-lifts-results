"""Solve and replay the first nonzero scalar degree using the original differential and SU(3) trace identities."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import json
from pathlib import Path

from su3_scalar_origin import GAUGINOS, add, formal, leading, read, reduce_image, scalar_count, serialize
from su3_subregular_image import write
import finite_trace_identity


def linear_combination(records, degree, denominator):
    if tuple(formal.ALPH) != tuple(finite_trace_identity.src.ALPH):
        raise ValueError("The letter tables disagree")
    gauginos = tuple(bytes((i,)) for i in sorted(GAUGINOS))
    bmn = {}
    identities = {}
    for kind, args, text in records:
        coefficient = denominator * Fraction(text)
        if not coefficient:
            raise ValueError("Zero coefficient in recipe")
        if kind == "classical":
            parent = tuple(bytes.fromhex(v) for v in args)
            add(bmn, parent, coefficient)
        elif kind == "trace_identity":
            first, second, rest = args
            spectators = tuple(bytes.fromhex(v) for v in rest)
            if any(formal.canon(t) != (1, t) for t in spectators):
                raise ValueError("Noncanonical spectator trace")
            pieces = (*gauginos, bytes.fromhex(first), bytes.fromhex(second))
            for v, sign in finite_trace_identity.identity(pieces).items():
                phase, term = formal.multi(v + spectators)
                if phase:
                    add(identities, term, coefficient*sign*phase)
        else:
            raise ValueError(f"Unknown recipe type {kind}")
    projected = dict(identities)
    for parent, coefficient in bmn.items():
        for t, n in formal.q(parent).items():
            if sum(map(scalar_count, t)) <= degree:
                add(projected, t, coefficient*n)
    if any(sum(map(scalar_count, t)) != degree for t in projected):
        raise ValueError("The proposed identities introduce a lower scalar layer")
    return bmn, identities, projected


def verify(level, directory, generate, next_residual):
    degree, inputs, denominator, total = leading(level, all_terms=True)
    recipes = json.loads((directory / f"j{level}_relations.json").read_bytes())
    bmn, identities, projected = linear_combination(recipes, degree, denominator)
    archive = directory / f"j{level}_primitive.json.gz"
    if generate:
        original = {w: c for w, c in total.items() if sum(map(scalar_count, w)) == degree}
        primitive, reduced_input = reduce_image(original)
        adjustment, reduced_relations = reduce_image(projected)
        if reduced_input != reduced_relations:
            raise ArithmeticError("The projected trace and classical relations do not give the input")
        for part in (bmn, adjustment):
            for w, c in part.items():
                add(primitive, w, -c)
    else:
        primitive = read(json.loads(gzip.decompress(archive.read_bytes())), denominator)

    for w, c in primitive.items():
        if len(b"".join(w)) != 10:
            raise ValueError("Incorrect parent length")
        for t, n in formal.q(w).items():
            add(total, t, c*n)
    for w, c in identities.items():
        add(total, w, -c)
    counts = Counter(sum(map(scalar_count, w)) for w in total)
    if any(n <= degree for n in counts):
        raise ArithmeticError("The full rational trace replay has a surviving low-order coefficient")
    if generate:
        write(archive, serialize(primitive, denominator))
    if next_residual:
        write(next_residual, serialize(total, denominator))
    dimensions = Counter(kind for kind, _, _ in recipes)
    return dict(level=level, input_by_scalar_degree=dict(sorted(inputs.items())), through_scalar_degree=degree,
                rational_primitive_terms=len(primitive), classical_relations=dimensions["classical"],
                four_slot_trace_identities=dimensions["trace_identity"], identity_polynomial_terms=len(identities),
                remaining_by_scalar_degree=dict(sorted(counts.items())), retained_residual=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, choices=(32, 34), required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--next-residual", type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("Enable assertions for the colour and cyclic sign checks")
    result = verify(args.level, args.data, args.generate, args.next_residual)
    out = args.data / f"j{args.level}_full_replay.json"
    if args.generate:
        out.write_text(json.dumps(result, indent=2)+"\n")
    elif json.loads(out.read_bytes()) != json.loads(json.dumps(result)):
        raise ArithmeticError("Replay and stored coordinates disagree")
    print(result)
