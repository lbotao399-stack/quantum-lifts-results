"""Evaluate the level-32 trace functional on every incident even parent trace."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

from itertools import product

from su3_level32_trace_factors import DATA, code, evaluator, multiply, trace
import block_first
import content_parents
import formal_trace_blocks as formal
from even_content_frames import conjugate
from finite_trace_identity import exterior_value


def incident(target):
    found=set()
    for i,splits in enumerate(formal.SPLITS):
        for _,pair in splits:
            left = Counter(target)
            left.subtract(pair)
            if min(left.values())<0:
                continue
            left[i]+=1
            found.add(bytes(sorted(left.elements())))
    return sorted(found)


def check(data, certificate):
    t=bytes(certificate["selected_content"])
    target=bytes(sorted(formal.ID[block_first.letter(k)] for k in t))
    translation=bytes([code(L) for L in formal.ALPH]+[0]*(256-len(formal.ALPH)))
    target_eval,point=evaluator(t,certificate["colour_subspace"],certificate["seed"])
    assert point["assignments"]==[ [k,i,tuple(m)] for k,i,m in certificate["assignments"] ]
    n=point["exterior_degree"]
    zero=((0,0,0),(0,0,0),(0,0,0))
    matrices={(k,i):zero for k in set(target) for i in (range(n) if formal.ODD[k] else [-1])}
    for k,i,m in point["assignments"]:
        matrices[formal.ID[block_first.letter(k)],i]=tuple(tuple(m[j:j+3]) for j in (0,3,6))
    basis = certificate["colour_subspace"]
    assert all(trace(v) == 0 for v in basis)
    gram = [[trace(multiply(x,y)) for y in basis] for x in basis]
    cartan = [trace(multiply(multiply(x,y),z))-trace(multiply(multiply(y,x),z))
              for x,y,z in product(basis,repeat=3)]
    assert gram == certificate["quadratic_gram"] and not any(map(any,gram))
    assert cartan == certificate["traces_of_commutators_times_basis"] and not any(cartan)
    value=Fraction()
    streamed=direct=0
    digest=hashlib.sha256()
    with gzip.open(data/"Q1.jsonl.gz","rb") as stream:
        for line in stream:
            digest.update(line)
            c,raw=json.loads(line)
            w=tuple(bytes(code(x) for x in tr) for tr in raw)
            if bytes(sorted(b"".join(w)))!=t:
                continue
            streamed+=1
            if any(len(tr)==2 for tr in w):
                continue
            native=tuple(bytes(formal.ID[block_first.letter(k)] for k in tr) for tr in w)
            coeff=exterior_value(native,matrices,n)
            assert coeff==target_eval(w)
            value+=Fraction(c)*coeff
            direct+=1
    assert digest.hexdigest()==certificate["source_sha256"] and direct==certificate["retained_words"]
    assert value==Fraction(certificate["Q1_source"]) and value
    def image(w):
        return sum(c*target_eval(tuple(tr.translate(translation) for tr in z))
                   for z,c in formal.q(w).items() if bytes(sorted(b"".join(z)))==target)
    by_shape=defaultdict(Counter)
    examples={}
    sizes=[]
    for g in incident(target):
        words=content_parents.all_words(g)
        sizes.append([list(g.translate(translation)),len(words)])
        for w in words:
            sign, partner=conjugate(w)
            if w>partner or (w==partner and sign<0):
                continue
            shape=tuple(sorted(map(len,w)))
            by_shape[shape]["tested"]+=1
            result=image(w)+sign*image(partner)
            if 2 in shape:
                assert result == 0
            if result:
                by_shape[shape]["nonzero"]+=1
                if shape not in examples:
                    examples[shape]=dict(word=[list(tr.translate(translation)) for tr in w],
                                         conjugate_word=[list(tr.translate(translation)) for tr in partner],
                                         conjugate_sign=sign, variation=result)
    return dict(Q1_source=str(value), input_content_words=streamed, independently_evaluated_words=direct,
                parent_contents=sizes, parents_by_trace_lengths=[dict(lengths=s,even_orbits=b["tested"],
                      nonzero=b["nonzero"],example=examples.get(s)) for s,b in sorted(by_shape.items())])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,default=DATA)
    p.add_argument("--certificate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not __debug__:
        p.error("verification requires enabled assertions")
    result=check(a.data,json.loads(a.certificate.read_bytes()))
    a.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result["Q1_source"], [(r["lengths"],r["nonzero"]) for r in result["parents_by_trace_lengths"]])
