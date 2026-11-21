"""Construct and check the level-32 first correction on the square-zero scalar locus."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from math import prod
from pathlib import Path

from su3_square_zero_image import ALPH, LETTER, ODD, SCALARS, adjacent_scalars, stream_rows
from fortuity.full.letters import Q_letter
from fortuity.full.words import canon_multi
from replay_formal_blocks import packed, expanded


DATA = Path(__file__).resolve().parent / "taylor_hooks/su3_level32"
GAUGINOS = {LETTER[1, 0, ()], LETTER[0, 1, ()]}
FERMIONS = {i for i, (p, m, t) in enumerate(ALPH) if p == m == 0 and len(t) == 2}
PROJECTOR = len(ALPH)
PARITY = ODD + (False,)


def canonical(trace):
    rotations = {}
    for i in range(len(trace)):
        left, right = trace[:i], trace[i:]
        sign = (-1)**(sum(PARITY[k] for k in left) * sum(PARITY[k] for k in right))
        rotated = right + left
        if rotations.setdefault(rotated, sign) != sign:
            return 0, b""
    first = min(rotations)
    return rotations[first], first


def multi(word):
    odd = [tr for tr in word if sum(PARITY[k] for k in tr) % 2]
    if len(set(odd)) != len(odd):
        return 0, ()
    return (-1)**sum(a > b for i, a in enumerate(odd) for b in odd[i+1:]), tuple(sorted(word))


def variation_parent(scalar, gaugino):
    p, m, _ = ALPH[gaugino]
    parent = p, m, ALPH[scalar][2]
    splits = {(u, v): n for n, u, v in Q_letter(parent)}
    sign = splits[ALPH[scalar], ALPH[gaugino]]
    if sign not in (-1, 1) or splits != {(ALPH[scalar], ALPH[gaugino]): sign,
                                         (ALPH[gaugino], ALPH[scalar]): -sign}:
        raise ValueError("unexpected classical derivative of a scalar")
    return LETTER[parent], sign


def commuting_gap(word):
    single = None
    for i, trace in enumerate(word):
        for start, letter in enumerate(trace):
            if letter not in SCALARS:
                continue
            left, right = trace[:start], trace[start:]
            rotated = right + left
            for j in range(1, len(rotated)):
                if rotated[j] in SCALARS:
                    rotation_sign = (-1)**(sum(ODD[k] for k in left) * sum(ODD[k] for k in right))
                    return i, rotated, j, rotation_sign, False
                if rotated[j] not in GAUGINOS:
                    break
            else:
                rotation_sign = (-1)**(sum(ODD[k] for k in left) * sum(ODD[k] for k in right))
                single = i, rotated, len(rotated), rotation_sign, True
    if single is None:
        raise ValueError("no scalar gap")
    return single


def terms():
    digest = hashlib.sha256()
    record = json.loads((DATA / "order1.json").read_bytes())
    result = Counter()
    count = 0
    for word, value in stream_rows(DATA / "Q1.jsonl.gz", digest):
        count += 1
        if not any(map(adjacent_scalars, word)):
            phase, canonical = canon_multi(expanded(word))
            result[packed(canonical)] += phase * Fraction(value)
    if digest.hexdigest() != record["output_sha256"] or count != record["direct_terms"]:
        raise ValueError("level-32 image disagrees with its archive")
    return {w: v for w, v in result.items() if v}


def build(source):
    correction = Counter()
    for word, value in source.items():
        letters = b"".join(word)
        if tuple(sum(k in group for k in letters) for group in (SCALARS, GAUGINOS, FERMIONS)) != (5, 2, 4):
            raise ValueError("unexpected surviving Taylor content")
        i, rotated, end, rotation_sign, single = commuting_gap(word)
        prefix_odd = sum(ODD[k] for trace in word[:i] for k in trace)
        for j in range(1, end):
            parent, split_sign = variation_parent(rotated[0], rotated[j])
            changed = rotated[1:j] + bytes([parent]) + rotated[j+1:]
            if single:
                changed = bytes([PROJECTOR]) + changed
            orientation = (-1)**(prefix_odd + j-1) * rotation_sign * split_sign
            inner_sign, inner = canonical(changed)
            if not inner_sign:
                continue
            phase, factors = multi(word[:i] + (inner,) + word[i+1:])
            if phase:
                correction[factors] -= value * orientation * inner_sign * phase
    return {w: v for w, v in correction.items() if v}


def save(path, correction):
    data = [[[tr.hex() for tr in w], str(v)] for w, v in sorted(correction.items())]
    raw = (json.dumps(data, separators=(",", ":"))+"\n").encode()
    path.write_bytes(gzip.compress(raw, mtime=0))


def reduce_trace(trace):
    if adjacent_scalars(trace):
        return 0, b""
    if PROJECTOR in trace:
        pos = trace.index(PROJECTOR)
        if trace[pos-1] in SCALARS:
            return 0, b""
        if trace[(pos+1) % len(trace)] in SCALARS:
            trace = trace[:pos] + trace[pos+1:]
    return canonical(trace)


def replay(source, path):
    data = json.loads(gzip.decompress(path.read_bytes()))
    variation = Counter()
    for w, value in data:
        word = tuple(bytes.fromhex(tr) for tr in w)
        prefix_sign = 1
        for i, trace in enumerate(word):
            for j, k in enumerate(trace):
                if k != PROJECTOR:
                    for n, u, v in Q_letter(ALPH[k]):
                        changed = trace[:j] + bytes([LETTER[u], LETTER[v]]) + trace[j+1:]
                        inner_sign, inner = reduce_trace(changed)
                        # Other factors may already contain the product of a projector and a scalar.
                        factors = [reduce_trace(t) for t in word[:i]] + [(inner_sign,inner)] + [reduce_trace(t) for t in word[i+1:]]
                        if all(s for s, _ in factors):
                            phase, target = multi(tuple(t for _, t in factors))
                            variation[target] += Fraction(value) * n * prefix_sign * phase * prod(s for s, _ in factors)
                if PARITY[k]:
                    prefix_sign = -prefix_sign
    residual = Counter(variation)
    residual.update(source)
    residual = {w: n for w, n in residual.items() if n}
    if residual:
        raise ValueError(f"{len(residual)} nonzero residual trace coefficients")
    return dict(level=32, equation="Q1 O0 + Q0 O1", full_first_image_terms=331201,
                terms_without_adjacent_scalars=len(source), correction_trace_terms=len(data),
                correction_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                correction_terms_with_projector=sum(any(PROJECTOR in bytes.fromhex(t) for t in w) for w, _ in data),
                Q0_terms_after_scalar_and_projector_relations=sum(bool(v) for v in variation.values()),
                residual_terms=len(residual))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--correction", type=Path, required=True)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", type=Path, required=True)
    args = parser.parse_args()
    source = terms()
    if args.build:
        save(args.correction, build(source))
    result = replay(source, args.correction)
    args.check.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)
