"""Extract the five four-spectator coefficients of the physical level-34 first image."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations
from pathlib import Path

import formal_trace_blocks as formal
import replay_formal_blocks as replay
from su3_levi_kernel import gzip_json


SCALAR = formal.ID[0, 0, (0,)]
EXTRA = (None, formal.ID[0, 0, (1,)], formal.ID[0, 0, (2,)],
         formal.ID[0, 1, ()], formal.ID[1, 0, ()])
NAMES = ("four_phi1", "three_phi1_phi2", "three_phi1_phi3",
         "three_phi1_minus", "three_phi1_plus")


@lru_cache(maxsize=150000)
def options(trace, extra):
    scalars = [k for k, letter in enumerate(trace) if letter == SCALAR]
    places = [(None, 1)]
    if extra is not None:
        places = [(k, (-1)**(formal.ODD[extra]*formal.par(trace[:k])))
                  for k, letter in enumerate(trace) if letter == extra]
    by_term = Counter()
    for place, sign in places:
        for n in range(min(4, len(scalars))+1):
            for slots in combinations(scalars, n):
                active = bytes(l for k, l in enumerate(trace) if k != place and k not in slots)
                if active:
                    phase, canonical = formal.canon(active)
                else:
                    phase, canonical = 2 + (-2)**len(trace), b""
                if phase:
                    by_term[n, canonical] += sign * phase
    return tuple((n, tr, c) for (n, tr), c in by_term.items() if c)


def coefficient(word, extra, out, scale):
    count = 4 if extra is None else 3
    placements = [(-1, 1)]
    if extra is not None:
        placements = []
        before = 0
        for i, trace in enumerate(word):
            if extra in trace:
                placements.append((i, (-1)**(formal.ODD[extra]*before)))
            before += formal.par(trace)
    for marked, outer_sign in placements:
        opts = [options(trace, extra if i == marked else None) for i, trace in enumerate(word)]
        def append(i, needed, traces, c):
            if i == len(opts):
                if not needed:
                    sign, canonical = formal.multi(traces)
                    if sign:
                        out[canonical] += c*sign
                return
            for n, trace, multiplicity in opts[i]:
                if n <= needed:
                    append(i+1, needed-n, traces + (trace,) if trace else traces, c*multiplicity)
        append(0, count, (), scale*outer_sign)


def extract(path):
    digest = hashlib.sha256()
    outputs = [Counter() for _ in EXTRA]
    count = 0
    for word, value in replay.stream_rows(path, digest):
        count += 1
        letters = b"".join(word)
        nscalar = letters.count(SCALAR)
        if nscalar < 3:
            continue
        scaled = -Fraction(value) * formal.DEN
        assert scaled.denominator == 1
        for extra, dest in zip(EXTRA, outputs):
            if (nscalar >= 4 if extra is None else extra in letters):
                coefficient(word, extra, dest, scaled.numerator)
    assert count == 517468 and digest.hexdigest() == formal.INPUT
    records = {}
    for name, dest in zip(NAMES, outputs):
        rows = [[[t.hex() for t in w], c] for w, c in sorted(dest.items()) if c]
        assert all(sum(len(t)//2 for t in w) == 7 for w, _ in rows)
        records[name] = rows
        print(name, len(rows), flush=True)
    return dict(source="Q1 O34", negative_input_sha256=formal.INPUT, input_terms=count,
                denominator=formal.DEN, coefficients=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(formal.__file__).with_name("negative_Q1.jsonl.gz"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    gzip_json(args.out, extract(args.input))
