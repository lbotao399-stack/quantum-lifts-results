"""Certify lowest-weight graviton products by all three BMN raising actions on a Cartan."""
import argparse
from collections import Counter
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path

from flint import nmod_mat

import mixed_incoming as m
from rank_two_bmn_product_replay import action, clean, highest, labelled, witness
from sp3_omega_product_replay import check_factor, odd_jet, odd_slot, ordered
from su3_bmn_parity_products import labelled as su3_labelled


@lru_cache(maxsize=24000)
def raise_word(slots, axis):
    result = Counter()
    sign = 1
    for i, slot in enumerate(slots):
        for j, jet in enumerate(slot):
            plus, minus, mask = jet
            if plus or minus or not mask:
                raise ValueError("Expected derivative-free BMN letters")
            if not mask & (1<<axis):
                target = (0,0,mask | (1<<axis))
                phase, word = ordered(slot[:j]+(target,)+slot[j+1:], odd_jet)
                if phase:
                    other = slots[:i]+(word,)+slots[i+1:]
                    result[other] += sign*phase*(-1)**((mask>>(axis+1)).bit_count())
            if odd_jet(jet):
                sign = -sign
    return tuple(clean(result).items())


def raise_polynomial(polynomial, axis):
    result = Counter()
    for w,n in polynomial.items():
        for v,c in raise_word(w,axis):
            result[v] += n*c
    return clean(result)


def top(polynomial):
    for axis in (0,1,2):
        polynomial = raise_polynomial(polynomial,axis)
    return polynomial


def check_action():
    labels = [(0,0,k) for k in range(1,8)]
    checked = 0
    for x in labels:
        for y in labels:
            row = {((x,),(y,)):1}
            for i in range(3):
                for j in range(i,3):
                    result = Counter(raise_polynomial(raise_polynomial(row,i),j))
                    result.update(raise_polynomial(raise_polynomial(row,j),i))
                    if clean(result):
                        raise ArithmeticError("BMN raising supercharges did not anticommute")
            output = top(row)
            for i,j in ((0,1),(1,0),(2,1),(1,2)):
                if action(output,i,j) != top(action(row,i,j)):
                    raise ArithmeticError("Triple raising did not intertwine the flavour roots")
            checked += 1
    return checked


def independent(rows):
    columns = {w:i for i,w in enumerate(sorted({w for row in rows for w in row}))}
    if not columns:
        return [],None
    table = [[0]*len(columns) for row in rows]
    for i,row in enumerate(rows):
        for w,n in row.items():
            table[i][columns[w]] = (n.numerator*pow(n.denominator,-1,m.P))%m.P
    keep = m.pivotcols(nmod_mat(table,m.P).transpose())
    return keep,witness([rows[i] for i in keep]) if keep else None


def compute(data, group, shape):
    unitary = group == "SU(3)"
    rows = []
    for record in data["rows" if unitary else "even_length_records"]:
        length = record["length"]
        products = [tuple(map(tuple,r)) for r in record["inputs" if unitary else "products"]]
        if not products:
            continue
        if not unitary:
            for spec in {factor for row in products for factor in row}:
                check_factor(spec)
        raised = []
        for row in products:
            coefficients = su3_labelled(row,True) if unitary else labelled(row)
            output = highest(top(coefficients),tuple(n+1 for n in shape))
            raised.append(output)
        keep, minor = independent(raised)
        rows.append(dict(length=length,trial_products=len(products),independent_input_indices=keep,
                         rank_lower=len(keep),minor=minor))
    return dict(group=group,weight=shape,records=rows,letter_pairs_checked=check_action())


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, nargs=3, required=True)
    parser.add_argument("--group",choices=("SU(3)","USp(4)"),required=True)
    parser.add_argument("--products",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:
        parser.error("Enable checks on the single traces and determinant")
    raw=args.products.read_bytes()
    data=json.loads(gzip.decompress(raw) if args.products.suffix==".gz" else raw)
    if isinstance(data,list):
        data,=[row for row in data if row["weight"]==args.weight]
    if data["weight"]!=args.weight or data["group"]!=args.group:
        raise ValueError("Product data have the wrong group or weight")
    result=compute(data,args.group,tuple(args.weight))
    result["input_sha256"]=hashlib.sha256(raw).hexdigest()
    args.out.write_text(json.dumps(result,separators=(",",":"))+"\n")
