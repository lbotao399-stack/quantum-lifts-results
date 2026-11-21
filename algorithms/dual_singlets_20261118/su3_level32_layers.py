"""Extend the two five-scalar parent solutions through the full six- and seven-scalar layers."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path

from su3_level32_content import DATA, formal, joint, physical
from su3_level32_symmetry import classical_controls, transform
from verify import add_scaled


def content(word):
    return bytes(sorted(b"".join(word)))


def first_targets():
    from g2_full_classical_screen_v1 import contents
    return sorted(bytes(sorted(formal.ID[letter] for letter in raw))
                  for raw in contents((1,1,7,3,3),11))


@lru_cache(100000)
def transform_trace(trace, kind):
    sign = 1
    letters = []
    for index in trace:
        factor, letter = transform(formal.ALPH[index],kind)
        letters.append(formal.ID[letter])
        sign *= factor
    cyclic, canonical = formal.canon(bytes(letters))
    assert cyclic
    return sign*cyclic, canonical


def transform_vector(vector, kind):
    result = {}
    for word, coefficient in vector.items():
        sign = 1
        traces = []
        for trace in word:
            factor, converted = transform_trace(trace,kind)
            sign *= factor
            traces.append(converted)
        reorder, canonical = formal.multi(traces)
        assert reorder and canonical not in result
        result[canonical] = sign*reorder*coefficient
    return result


def images(vector):
    derivatives = transform_vector(vector,"derivatives")
    flavours = transform_vector(vector,"flavours")
    both = transform_vector(derivatives,"flavours")
    assert both == transform_vector(flavours,"derivatives")
    return [vector,derivatives,flavours,both]


def read_correction(path):
    saved = joint.s.read(path)
    assert saved["rank"]["kind"] == "exact_solution" and saved["parent_scalar_count"] == 5
    result = {joint.s.unpack_hex(w):Fraction(c) for w,c in saved["correction"]}
    assert all(content(w).count(0) == 5 for w in result)
    return result


def assemble(pair, four):
    result = {}
    for path, divisor in ((pair,4),(four,1)):
        for row, sign in zip(images(read_correction(path)),(1,-1,-1,1),strict=True):
            add_scaled(result,row,Fraction(sign,divisor))
    for kind in ("derivatives","flavours"):
        assert transform_vector(result,kind) == {w:-c for w,c in result.items()}
    return result


def residual(data, correction, keep):
    digest = hashlib.sha256()
    rhs = defaultdict(dict)
    targets = set(first_targets())
    for k in range(len(formal.SPLITS)):
        initial = int(k==0)
        assert all(initial <= p.count(0) <= initial+1 for _,p in formal.SPLITS[k])
    assert not formal.SPLITS[0]
    with gzip.open(data/"Q1.jsonl.gz","rb") as stream:
        for line in stream:
            digest.update(line)
            c, raw = json.loads(line)
            word = tuple(bytes(formal.ID[a,b,tuple(f)] for a,b,f in tr) for tr in raw)
            t = content(word)
            assert t in targets
            if keep(t):
                rhs[t][word] = Fraction(c)
    info = json.loads((data/"order1.json").read_bytes())
    assert digest.hexdigest() == info["output_sha256"]
    for word, coefficient in correction.items():
        for v,n in formal.q(word).items():
            t = content(v)
            assert t in targets
            if keep(t):
                add_scaled(rhs[t],{v:n},coefficient)
    return {t:rhs[t] for t in sorted(targets) if keep(t)}


def certify(data, correction, engine, work):
    work.mkdir(exist_ok=True,parents=True)
    rows = residual(data,correction,lambda t:t.count(0)>=6)
    checks = {}
    for t,r in rows.items():
        if r:
            check = physical.polynomials([r],engine,work)[0]
            assert not check["nonzero_monomials"], (t.hex(),check)
            checks[t.hex()] = check
    counts = Counter(t.count(0) for t in rows)
    assert counts == {6:78,7:13}
    return dict(target_contents_by_scalar_count=dict(sorted(counts.items())),formal_zero=len(rows)-len(checks),
                polynomial_checks=checks, total_collected_monomials=sum(v["collected_monomials"] for v in checks.values()),
                total_nonzero_monomials=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,default=DATA)
    parser.add_argument("--pair",required=True,type=Path)
    parser.add_argument("--four",required=True,type=Path)
    parser.add_argument("--engine",required=True,type=Path)
    parser.add_argument("--work",required=True,type=Path)
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    args.work.mkdir(exist_ok=True,parents=True)
    for kind in ("derivatives","flavours"):
        assert classical_controls(kind) == 180
    correction = assemble(args.pair,args.four)
    result = certify(args.data,correction,args.engine.resolve(),args.work)
    result.update(correction_terms=len(correction),
                  input_sha256={k:hashlib.sha256(p.read_bytes()).hexdigest()
                                for k,p in (("pair",args.pair),("four",args.four))})
    joint.s.write(args.work/"correction.json.gz",[[joint.s.hexes(w),str(c)] for w,c in sorted(correction.items())])
    joint.s.write(args.work/"six_and_seven.json",result)
    print(len(correction),result["formal_zero"],result["total_collected_monomials"],flush=True)


if __name__=="__main__":
    main()
