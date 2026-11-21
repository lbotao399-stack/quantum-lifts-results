"""Project the level-32 first correction under derivative and equal-weight flavour exchange."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path

from su3_level32_walls import DATA, algebra, block_first
from fortuity.full.letters import Q_letter, letters_up_to
from verify import add_scaled


def transform(letter, kind):
    plus,minus,flavours = letter
    if kind == "derivatives":
        return 1, (minus,plus,flavours)
    changed = tuple((0,2,1)[i] for i in flavours)
    inversions = sum(i>j for k,i in enumerate(changed) for j in changed[k+1:])
    return (-1)**inversions,(plus,minus,tuple(sorted(changed)))


@lru_cache(80000)
def transform_trace(tr, kind):
    factor = 1
    result = []
    for code in tr:
        sign, letter = transform(block_first.letter(code),kind)
        result.append(letter)
        factor *= sign
    sign, trace = algebra.fw.canon_word(tuple(result))
    assert sign
    return factor*sign,trace


def transform_word(word,kind):
    factor = 1
    traces = []
    for tr in word:
        sign, result = transform_trace(tr,kind)
        factor *= sign
        traces.append(result)
    sign,key = algebra.fw.canon_multi(tuple(traces))
    assert sign
    return factor*sign,tuple(bytes(map(block_first.code,tr)) for tr in key)


def classical_controls(kind):
    count = 0
    for letter in letters_up_to(1,1):
        sign, changed = transform(letter,kind)
        direct = Counter()
        for c,left,right in Q_letter(changed):
            direct[left,right] += sign*c
        other = Counter()
        for c,left,right in Q_letter(letter):
            s, l = transform(left,kind)
            t, r = transform(right,kind)
            other[l,r] += s*t*c
            count += 1
        assert direct == other
    return count


def project(data, certificate):
    info = json.loads((data/"order1.json").read_bytes())
    digest = hashlib.sha256()
    rows = {}
    with gzip.open(data/"Q1.jsonl.gz","rb") as stream:
        for line in stream:
            digest.update(line)
            c,raw=json.loads(line)
            word=tuple(bytes(map(block_first.code,tr)) for tr in raw)
            value=Fraction(c)*info["common_denominator"]
            assert value.denominator == 1 and word not in rows
            rows[word]=value.numerator
    assert digest.hexdigest() == info["output_sha256"] and len(rows) == info["direct_terms"]
    result={tuple(map(bytes,w)):Fraction(c) for w,c in certificate["correction"]}
    controls={}
    for kind in ("derivatives","flavours"):
        splits=classical_controls(kind)
        sample,first=next(iter(rows.items()))
        sign,changed=transform_word(sample,kind)
        eigen=Fraction(rows[changed],first*sign)
        assert eigen in (-1,1)
        for w,c in rows.items():
            sign,changed=transform_word(w,kind)
            assert rows[changed] == int(eigen)*sign*c
        other={}
        for w,c in result.items():
            sign,changed=transform_word(w,kind)
            other[changed]=sign*c
        add_scaled(result,other,eigen)
        result={w:c/2 for w,c in result.items()}
        controls[kind]=dict(eigenvalue=int(eigen),classical_split_terms=splits,checked_quantum_terms=len(rows))
    for kind, control in controls.items():
        for w,c in result.items():
            s,k=transform_word(w,kind)
            assert result[k] == s*control["eigenvalue"]*c
    return dict(kind="solution",first_image_sha256=digest.hexdigest(),original_terms=len(certificate["correction"]),
                exchanges=controls,correction=[[list(map(list,w)),str(c)] for w,c in sorted(result.items())])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--certificate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    result=project(a.data,json.loads(a.certificate.read_bytes()))
    a.out.write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(result["exchanges"],len(result["correction"]),flush=True)
