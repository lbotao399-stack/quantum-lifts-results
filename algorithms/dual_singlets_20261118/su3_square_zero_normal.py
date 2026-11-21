"""Test the level-32 invariant correction modulo the square of the scalar-product ideal."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from math import comb
from pathlib import Path

from su3_square_zero_primitive import ALPH, DATA, SCALARS, stream_rows
import formal_trace_blocks as formal


@lru_cache(maxsize=30000)
def disjoint_pairs(trace):
    """Maximum number of disjoint cyclically adjacent pairs of scalar letters."""
    cut = next((i for i,k in enumerate(trace) if k not in SCALARS), None)
    if cut is None:
        return len(trace)//2
    count = run = 0
    for k in trace[cut+1:] + trace[:cut+1]:
        if k in SCALARS:
            run += 1
        else:
            count += run//2
            run = 0
    return count


def canonical(word):
    traces = []
    sign = 1
    for tr in word:
        phase, norm = formal.canon(tr)
        traces.append(norm)
        sign *= phase
    reorder,norm = formal.multi(traces)
    return sign*reorder,norm


def filtered_residual(correction):
    assert tuple(formal.ALPH) == ALPH
    digest = hashlib.sha256()
    result = defaultdict(Counter)
    totals = Counter()
    def insert(word,c):
        t = bytes(sorted(b"".join(word)))
        totals[t] += 1
        if sum(map(disjoint_pairs,word))<2:
            sign,w = canonical(word)
            if sign:
                result[t][w] += sign*Fraction(c)
    count = 0
    for w,c in stream_rows(DATA/"Q1.jsonl.gz",digest):
        insert(w,c)
        count += 1
    record = json.loads((DATA/"order1.json").read_bytes())
    assert count==record["direct_terms"] and digest.hexdigest()==record["output_sha256"]
    data = json.loads(gzip.decompress(correction.read_bytes()))
    for encoded,c in data:
        w = tuple(bytes.fromhex(t) for t in encoded)
        for tr,n in formal.q(w).items():
            insert(tr,n*Fraction(c))
    return {t:{w:c for w,c in result[t].items() if c} for t in totals}


def unit(i,j):
    return tuple(int(k==3*i+j) for k in range(9))


def add(a,b,scale):
    return tuple(x+scale*y for x,y in zip(a,b,strict=True))


def transverse_value(t,row,seed):
    from su3_level32_trace_factors import code, evaluator

    trans = bytes(map(code,ALPH)) + bytes(256-len(ALPH))
    colour_basis = (add(unit(0,0),unit(2,2),-1),add(unit(1,1),unit(2,2),-1),
                    *(unit(i,j) for i in range(3) for j in range(3) if i!=j))
    base = (unit(0,1), unit(0,2), add(unit(0,1),unit(0,2),1))
    direction = (unit(1,0),unit(2,1),add(unit(0,0),unit(1,1),-1))
    power = sum(k in SCALARS for k in t)
    values = []
    assignment = None
    for v in range(power+1):
        fixed = {trans[k]:add(base[ALPH[k][2][0]],direction[ALPH[k][2][0]],v) for k in t if k in SCALARS}
        ev,info = evaluator(t.translate(trans),colour_basis,seed,fixed)
        if v==0:
            assignment = info
        values.append(sum(c*ev(tuple(tr.translate(trans) for tr in w)) for w,c in row.items()))
    assert values[0]==0
    derivative = sum(Fraction((-1)**(k+1)*comb(power,k),k)*(values[k]-values[0]) for k in range(1,power+1))
    return dict(content=t.hex(),input_letters=[ALPH[k] for k in t],seed=seed,degree_bound=power,
                retained_trace_terms=len(row),values=list(map(str,values)),linear_coefficient=str(derivative),
                scalar_base=base,scalar_direction=direction,assignment_at_zero=assignment)


def run(correction):
    rows = filtered_residual(correction)
    sizes = defaultdict(Counter)
    for t,z in rows.items():
        counts = sizes[sum(k in SCALARS for k in t)]
        counts["contents"] += 1
        counts["remaining_trace_terms"] += len(z)
        counts["zero_by_disjoint_pairs"] += not z
    for t in sorted(rows,key=lambda k:(len(rows[k]), k)):
        if not rows[t]:
            continue
        witness = transverse_value(t,rows[t],146087)
        if Fraction(witness["linear_coefficient"]):
            return dict(level=32,correction_sha256=hashlib.sha256(correction.read_bytes()).hexdigest(),
                        residual_by_scalar_degree=sizes,contents=len(rows),transverse_witness=witness)
    raise ValueError("this family of transverse directions gave no nonzero coefficient")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--correction",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a = p.parse_args()
    if not __debug__:p.error("enable assertions")
    a.out.write_text(json.dumps(run(a.correction),indent=2)+"\n")
