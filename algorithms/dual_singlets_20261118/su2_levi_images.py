"""Export the published SU(3) first variations in the native two-colour Taylor convention."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
from io import StringIO
import json
from math import prod
from pathlib import Path

import formal_trace_blocks as formal


HERE = Path(__file__).resolve().parent
INPUTS = {
    32: (HERE/"taylor_hooks/su3_level32/Q1.jsonl.gz",
         "f7810267d100f7cc9c9b4c93e045540ab2291df8b8bb35eeb4ad77feddd47a8a", 3, 331201),
    34: (HERE.parent/"su3_j34_20261118/negative_Q1.jsonl.gz",
         formal.INPUT, 120, 517468),
}


def native(letter):
    r, s, flavours = letter
    if len(flavours) == 0:
        return (3, 0, r, s), -1
    if len(flavours) == 1:
        return (0, flavours[0], r, s), 1
    if len(flavours) == 2:
        missing = 3 - sum(flavours)
        return (1, missing, r, s), -1 if missing == 1 else 1
    return (2, 0, r, s), 1


def check_convention():
    """Expand each native supercommutator and compare ordered products with every input split."""
    for k, letter in enumerate(formal.ALPH):
        (kind, f, r, s), phase = native(letter)
        left = Counter()
        for sign, (i, j) in formal.SPLITS[k]:
            x, a = native(formal.ALPH[i])
            y, b = native(formal.ALPH[j])
            left[x, y] += sign*a*b
        if kind == 3:
            brackets = [(3, 0, 3, 0, Fraction(-1, 2))]
        else:
            brackets = [(3, 0, kind, f, -1)]
            if kind == 1:
                brackets.append((0, (f+1)%3, 0, (f+2)%3, 1))
            if kind == 2:
                brackets.extend((0, i, 1, i, -1) for i in range(3))
        right = Counter()
        for x, i, y, j, c in brackets:
            for u in range(r+1):
                for v in range(s+1):
                    if (x == 3 and u+v == 0) or (y == 3 and r+s-u-v == 0):
                        continue
                    first, second = (x, i, u, v), (y, j, r-u, s-v)
                    right[first, second] += phase*c
                    right[second, first] += phase*c*(1 if x in (1,3) and y in (1,3) else -1)
        assert left == right


def export(level, output):
    path, expected, denominator, count = INPUTS[level]
    trace_ids = {}
    traces = []
    signs = []
    digest = hashlib.sha256()
    terms = StringIO()
    total = 0
    alphabet = [native(l) for l in formal.ALPH]
    with gzip.open(path, "rb") as stream:
        for line in stream:
            digest.update(line)
            c, row = json.loads(line)
            ids = []
            for trace in row:
                word = bytes(formal.ID[r, s, tuple(f)] for r, s, f in trace)
                if word not in trace_ids:
                    trace_ids[word] = len(traces)
                    signs.append(prod(alphabet[i][1] for i in word))
                    traces.append(word)
                ids.append(trace_ids[word])
            # The level-34 stream stores minus Q1. Both exports are plus Q1.
            value = Fraction(c)*denominator*prod(signs[i] for i in ids)*(1 if level == 32 else -1)
            assert value.denominator == 1 and abs(value) < 2**63
            print(int(value), len(ids), *ids, file=terms)
            total += 1
    assert digest.hexdigest() == expected and total == count
    with output.open("w") as out:
        print(len(alphabet), len(traces), total, file=out)
        for coords, _ in alphabet:
            print(*coords, file=out)
        for word in traces:
            print(len(word), *word, file=out)
        out.write(terms.getvalue())
    return dict(level=level, source_terms=total, traces=len(traces), denominator=denominator,
                input_sha256=expected, exported_sha256=hashlib.sha256(output.read_bytes()).hexdigest())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level", type=int, choices=INPUTS, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if not __debug__:
        p.error("convention and input checks require enabled assertions")
    check_convention()
    print(json.dumps(export(a.level, a.out), sort_keys=True))
