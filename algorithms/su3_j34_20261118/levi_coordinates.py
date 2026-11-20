"""Test the three coupled scalar coefficients on genuine SU(3) highest weights.

The diagonal scalar is diag(1+x,-1+x,-2x). The other five fields
occupy the recorded diagonal, upper, and lower entries in the first two
colours. Direct three-colour multiplication recovers the even powers
0, 2, 4 by interpolation of a polynomial of degree at most five.
"""
import argparse
import hashlib
import json
import random
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from sympy import Matrix
import formal_trace_blocks as f
import replay_formal_blocks as replay
from block_first import code, letter
from g2_full_classical_screen_v1 import sector
import su_trace_adapter_v1 as su


CONTENT = tuple(sorted([1]*5 + [3, 5, 6, 6, 8, 18]))
TARGET = bytes(sorted(f.ID[letter(k)] for k in CONTENT))
PARENT_CHARGES = sector((1, 1, 7, 4, 3), 10)
COLOUR_CODES = tuple(code(l) for l in f.ALPH)
VALUES = (0, 1, -1, 2, -2)


@lru_cache(65536)
def trace_at(trace):
    # Encodings 63 and 64 distinguish the two odd entries of letter 6.
    out = [0]*len(VALUES)
    for start in range(3):
        pos = start
        signs = 1
        powers = [0, 0, 0]
        for k in trace:
            if k == 1:
                powers[pos] += 1
            elif k in (5, 18):
                signs *= (1, -1, 0)[pos]
            elif k in (8, 63) and pos == 0:
                pos = 1
            elif k in (3, 64) and pos == 1:
                pos = 0
            else:
                signs = 0
                break
        if signs and pos == start:
            a,b,c = powers
            for j,x in enumerate(VALUES):
                out[j] += signs*(1+x)**a*(-1+x)**b*(-2*x)**c
    return tuple(out)


@lru_cache(60000)
def coordinates(word):
    sums = [0]*len(VALUES)
    for order in ((63, 64), (64, 63)):
        chosen = iter(order)
        exterior = []
        product = [1]*len(VALUES)
        slots = {3: 11, 5: 15, 63: 19, 64: 20, 8: 25}
        for tr in word:
            row = tuple(next(chosen) if k == 6 else k for k in map(COLOUR_CODES.__getitem__, tr))
            exterior.extend(slots[k] for k in row if k in slots)
            product = [a*b for a,b in zip(product, trace_at(row))]
        assert sorted(exterior) == [11, 15, 19, 20, 25]
        sign = (-1)**sum(i>j for p,i in enumerate(exterior) for j in exterior[p+1:])
        sums = [a+sign*b for a,b in zip(sums, product)]
    c0, pos, neg, two, minus = sums
    assert (pos+neg)%2 == (two+minus)%2 == 0
    unit, doub = (pos+neg)//2-c0, (two+minus)//2-c0
    assert (doub-4*unit)%12 == 0
    c4 = (doub-4*unit)//12
    return c0, unit-c4, c4


def incident():
    target = Counter(TARGET)
    return {bytes(sorted(list((target-Counter(pair)).elements())+[i]))
            for i, terms in enumerate(f.SPLITS) for _,pair in terms if Counter(pair) <= target}


@lru_cache(32768)
def image(w):
    out = [0, 0, 0]
    for z, c in f.q(w).items():
        if bytes(sorted(b"".join(z))) == TARGET:
            out = [a+c*b for a,b in zip(out, coordinates(z))]
    return out


def highest(w, parents):
    expr = su.projection.project(su.parity_project({replay.expanded(w): 1}, 1), PARENT_CHARGES)
    for root in ((0,1), (1,2), ("spin",1)):
        assert not su.projection.action(expr, *root)
    assert su.conjugate(expr) == expr
    out = [Fraction(0)]*3
    for v, c in expr.items():
        u = replay.packed(v)
        if bytes(sorted(b"".join(u))) in parents:
            out = [a+c*b for a,b in zip(out, image(u))]
    return len(expr), out


def search(parents):
    rng = random.Random(342011)
    shapes = [(2,2,2,2,2), (2,2,3,3), (3,3,4), (4,6), (5,5), (10,)]
    selected, rows = [], []
    for attempt in range(4000):
        letters = list(rng.choice(sorted(parents)))
        rng.shuffle(letters)
        traces = []
        for n in shapes[attempt%len(shapes)]:
            sg, tr = f.canon(bytes(letters[:n]))
            if not sg:
                break
            traces.append(tr)
            del letters[:n]
        else:
            sg, w = f.multi(tuple(traces))
            if not sg or not any(image(w)):
                continue
            count, v = highest(w, parents)
            if Matrix(rows+[v]).rank() > len(rows):
                selected.append((w, count))
                rows.append(v)
                if len(rows) == 3:
                    return selected, rows
    raise ValueError("Search did not span the three coordinates")


def run(rhs, functional, recipes):
    saved = json.loads(functional.read_bytes())
    terms = saved["functional"]
    assert [e[0][:2] for e in terms] == [[0,0],[1,2],[1,4]]
    for (s, p, content, diagonal, key), _ in terms:
        assert tuple(sorted(content+[1]*p)) == CONTENT and diagonal == 1
        assert int(key[0]) == 1 << (4*54) and int(key[1]) == sum(1 << i for i in (11,15,19,20,25))
    functional_row = [Fraction(c) for _,c in terms]
    target = [Fraction(0)]*3
    h = hashlib.sha256()
    rhs_terms = 0
    for w, c in replay.stream_rows(rhs, h):
        if bytes(sorted(b"".join(w))) == TARGET:
            rhs_terms += 1
            target = [a+Fraction(c)*b for a,b in zip(target, coordinates(w))]
    assert h.hexdigest() == saved["input"] == f.INPUT
    pairing = sum(a*b for a,b in zip(functional_row,target))
    assert pairing == Fraction(saved["pairing"]) == Fraction(1,120)
    parents = incident()
    if recipes:
        words = [replay.unpack_hex(r["word"]) for r in json.loads(recipes.read_bytes())["recipes"]]
        tests = [highest(w, parents) for w in words]
        selected = list(zip(words, [n for n,_ in tests]))
        rows = [v for _,v in tests]
    else:
        selected, rows = search(parents)
    mat = Matrix(rows).T
    assert mat.det()
    coefficients = mat.inv()*Matrix(target)
    assert mat*coefficients == Matrix(target)
    return dict(rhs_sha256=h.hexdigest(), target_content=TARGET.hex(),
                central_powers=[0,2,4], three_colour_rhs_terms=rhs_terms,
                rhs=list(map(str,target)), paired_rhs=str(pairing), incident_parent_contents=len(parents),
                highest_weight_doubled_charges=PARENT_CHARGES,
                recipes=[dict(word=replay.hexes(w), projected_words=n, coordinates=list(map(str,v)),
                              coefficient=str(c), functional_pairing=str(sum(a*b for a,b in zip(v,functional_row))))
                         for (w,n),v,c in zip(selected,rows,coefficients)],
                determinant=str(mat.det()), residual=["0"]*3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rhs", type=Path, default=f.HERE/"negative_Q1.jsonl.gz")
    parser.add_argument("--functional", type=Path, default=f.HERE/"levi"/"joint.json")
    parser.add_argument("--recipes", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    assert __debug__
    doc = run(args.rhs, args.functional, args.recipes)
    args.out.write_text(json.dumps(doc, indent=2)+"\n")
    print(doc["rhs"], doc["determinant"], doc["residual"])
