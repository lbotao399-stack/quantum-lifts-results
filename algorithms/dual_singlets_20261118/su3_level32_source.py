"""Expand the level-32 SU(3) operator from the published arXiv:2412.08695v2 TeX."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
from itertools import permutations, product
import json
from pathlib import Path
import re
import sys


def equation(paper):
    # The last occurrence of this label is the printed v2 formula. The earlier
    # occurrence is inside an inactive TeX conditional.
    text = paper.read_text().rsplit(r"\ba\label{32cohomology}", 1)[1].split(r"\ea", 1)[0]
    text = text.split("&=", 1)[1].replace(r"\nn", "").replace(r"\\", "").replace(r"\ .", "")
    text = "".join(text.split()).replace("&", "").replace("|", "")
    text = re.sub(r"\(r_\{(\d+)\}\^\{\(0,\d\)\}\)", r"r\1", text)
    text = text.replace("(", "").replace(")", "").replace(r"\ep", "E")
    text = re.sub(r"\{v\^\{([^}]*)\}\}", r"v^{\1}", text)
    terms = []
    for raw in re.finditer(r"([+-]\d+)([^+-]+)", text):
        coefficient, body = raw.groups()
        factors, balance = [], Counter()
        position = 0
        for match in re.finditer(r"(E|u|v|r\d+)(?:\^\{([^}]+)\})?(?:_\{([^}]+)\})?", body):
            assert match.start() == position
            position = match.end()
            symbol, upper, lower = match.groups()
            indices = []
            for group, sign in ((upper, 1), (lower, -1)):
                codes = re.findall(r"a_\d+|[bijkl]", group or "")
                assert "".join(codes) == (group or "")
                for code in codes:
                    balance[code] += sign
                indices.append(codes)
            up, down = indices
            if symbol in ("u", "v"):
                symbol += str(len(up) + len(down))
            elif symbol.startswith("r"):
                symbol += f"_{len(up)}{len(down)}"
            factors.append((symbol, " ".join(up+down)))
        assert position == len(body) and {name:n for name,n in balance.items() if n} == dict.fromkeys("ijkl", 1)
        terms.append((int(coefficient), factors))
    assert len(terms) == 47 and "".join(m.group() for m in re.finditer(r"([+-]\d+)([^+-]+)", text)) == text
    return terms


def additional_preimages(dmk):
    # App. A.1. Only the fundamental [1,0] and [0,2] tensors are missing
    # from the level-30 implementation. Their two lower labels are symmetrized.
    def tr(*letters):
        return dmk.tr(*({i: 1} for i in letters))

    @lru_cache(None)
    def vector(i):
        terms = []
        for a, b in product(range(3), repeat=2):
            A, B = 3+a, 3+b
            terms.extend([
                (3, tr(i,A,a,b,B)),
                (1, dmk.mul(tr(i,a),tr(b,A,B))), (1, dmk.mul(tr(i,a),tr(b,B,A))),
                (-3, dmk.mul(tr(i,A),tr(a,b,B))), (3, dmk.mul(tr(i,A),tr(b,a,B))),
                (-1, dmk.mul(tr(i,A,B),tr(a,b)))])
        return dmk.add(*terms)

    @lru_cache(None)
    def tensor(i, j):
        terms = []
        for left, right in ((i,j), (j,i)):
            J = 3+right
            for a,b,c in permutations(range(3)):
                e = Fraction(5*dmk.eps(a,b,c), 18)
                terms.extend([(2*e, tr(3+left,J,a,b,c)), (e, tr(3+left,a,J,b,c))])
            for a,b in permutations(range(3),2):
                e = Fraction(dmk.eps(a,b,left), 2)
                if not e:
                    continue
                for c in range(3):
                    C = 3+c
                    for scalars in ((a,b,c), (a,c,b), (c,a,b)):
                        for k in range(4):
                            terms.append((e, tr(J,*scalars[:k],C,*scalars[k:])))
                    terms.extend([
                        (-5*e/3, dmk.mul(tr(J,a,b),tr(C,c))),
                        (-e/3, dmk.mul(tr(J,a,c),tr(C,b))),
                        (-e/3, dmk.mul(tr(J,c,a),tr(C,b))),
                        (e/3, dmk.mul(tr(J,b),tr(C,a,c))),
                        (e/3, dmk.mul(tr(J,b),tr(C,c,a)))])
        return dmk.add(*terms)
    return {"r14_10":vector, "r14_02":tensor}


def build(dmk, terms):
    dmk._RAW.update(additional_preimages(dmk))
    total = Counter()
    for k, (coefficient, factors) in enumerate(terms, 1):
        # This printed operator includes the signs of its odd prefactors.
        # Applying Q_0 restores those signs in the relation among the R tensors.
        expression = dmk.term_expr(factors, dict.fromkeys("ijkl", 0))
        for word, number in expression.items():
            total[word] += coefficient * number
        print(f"Expression {k}/47, {len(expression)} terms", flush=True)
    return {w:c for w,c in total.items() if c}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--fortuity-code", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    modules = args.fortuity_code.resolve()
    pins = {args.paper:"c9e328f5aeb008480507640d70b1a499ede31939f5889bb9513f9c341b764c4a",
            modules/"fortuity/bmn/dmk.py":"b955ec8feb8763fa111147a09fa6da02aed73188c0aa14d08fee46d400150823",
            modules/"fortuity/bmn/words.py":"e2ecc349c6485c0d2063d60ba955d820d7f5ba01032a84645caba66045de5593"}
    assert all(hashlib.sha256(path.read_bytes()).hexdigest()==value for path,value in pins.items())
    sys.path.insert(0, str(modules))
    from fortuity.bmn import dmk, words
    from su3_bmn_descendants import even_part, flavour_root
    source = build(dmk, equation(args.paper))
    assert words.expr_sector(source) == (5,8,4,4)
    assert not flavour_root(source, 0, 1, words) and not flavour_root(source, 1, 2, words)
    # The colour projection commutes with Q_0. The quartic scalar restriction
    # below detects the even component of this published class.
    source = even_part(source, words)
    record = dict(level=32, weight=[8,4,4], colour_parity="even", length=11, n_terms=len(source),
                  terms=[[str(c), w] for w,c in sorted(source.items())])
    data = json.dumps(record, separators=(",", ":")).encode() + b"\n"
    args.out.write_bytes(gzip.compress(data, mtime=0))
    print("Highest-weight even component", len(source), hashlib.sha256(data).hexdigest(), flush=True)
