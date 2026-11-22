"""Transport the explicit BMN first correction to the Cartan-null nine-letter target."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, product as pairs
import json
from pathlib import Path

import product
from first_correction import construct
from fortuity.full.letters import is_ferm
from q1d_exact import pair_terms


def pair_check(x,y,spin,flavour):
    terms = Counter()
    def change(a):
        return product.substitution(a,spin,flavour,True)
    for c,l,r in pair_terms(x,y):
        a,z = change(l)
        if a:
            terms[z,r] += a*c
        a,z = change(r)
        if a:
            terms[l,z] += a*c*(-1)**is_ferm(l)
    contributions = any(terms.values())
    a,z = change(x)
    if a:
        for c,l,r in pair_terms(z,y):
            terms[l,r] += a*c
    b,z = change(y)
    if b:
        for c,l,r in pair_terms(x,z):
            terms[l,r] += b*c*(-1)**is_ferm(x)
    assert not any(terms.values()), (x,y,spin,flavour)
    return contributions


def streamed_digest(polynomial):
    result = sha256(b"[")
    for i,w in enumerate(sorted(polynomial)):
        if i:
            result.update(b",")
        result.update(json.dumps((w,polynomial[w]),separators=(",",":")).encode())
    result.update(b"]")
    return result.hexdigest()


def verify():
    source,correction,record = construct()
    source,den = product.integral(source)
    assert den==record["correction_denominator"]==180
    assert streamed_digest(correction)==record["correction_numerator_sha256"]
    alphabet = [(a,b,t) for n in range(1,4) for a in range(n+1)
                for b in range(n+1-a) for t in combinations(range(3),n-a-b)]
    assert len(alphabet)==37
    stages = []
    for spin,flavour in product.STEPS:
        actual = {letter for w in source for trace in w for letter in trace}
        assert actual<=set(alphabet)
        classical = {letter for w in correction for trace in w for letter in trace}
        for letter in actual | classical:
            product.check_letter(letter,spin,flavour,True)
        active = sum(pair_check(x,y,spin,flavour) for x,y in pairs(alphabet,repeat=2))
        source = product.action(source,spin,flavour,True)
        correction = product.action(correction,spin,flavour,True)
        stages.append(dict(spin=spin,flavour=flavour,letters_checked=len(actual|classical),
                           ordered_pairs=len(alphabet)**2,nonzero_output_actions=active,
                           source_terms=len(source),correction_terms=len(correction)))
    source_hash = streamed_digest(source)
    data = json.loads(Path(__file__).with_name("product.json").read_bytes())
    assert source_hash==data["raised_numerator_sha256"]
    assert [s["source_terms"] for s in stages]==data["terms_after_substitutions"][:4]
    assert product.charge(correction)==((3,3,4,1,1),8)
    return dict(stages=stages,source_numerator_sha256=source_hash,common_denominator=den,
                correction_numerator_sha256=streamed_digest(correction))


if __name__=="__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",type=Path)
    args = p.parse_args()
    if not __debug__:
        p.error("enable assertions")
    encoded = json.dumps(verify(),indent=2)+"\n"
    if args.check:
        assert encoded==args.check.read_text()
    else:
        print(encoded,end="")
