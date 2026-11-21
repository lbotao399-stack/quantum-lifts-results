"""Compute the Cartan image of quadratic, cubic or quartic protected trace products in a Taylor weight."""
import argparse
import hashlib
from itertools import islice
import json
from pathlib import Path

import numpy as np
from flint import nmod_mat

from sp3_omega_product_replay import check_factor
from sp3_omega_products import all_products, cartan
from su3_bmn_parity_products import labelled
from taylor_hook_cohomology import weights


PRIME = 2147483647
SECTORS = {"su30":("SU3",3,9), "su32":("SU3",4,10), "sp36":("Sp2",6,12)}


def pivots(matrix):
    reduced, rank = matrix.rref(inplace=True)
    result = []
    j = 0
    for i in range(rank):
        while not reduced[i,j]:
            j += 1
        result.append(j)
        j += 1
    return result


def coordinates(product, group):
    return labelled(product,True) if group == "SU3" else cartan(product,2)


def matrix(rows, width, prime):
    a = np.zeros((len(rows),width),dtype=np.int64)
    for i,row in enumerate(rows):
        for k,v in row.items():
            a[i,k] = int(v)%prime
    return a


def record(case,index,raised_spin=False):
    group, tail, length = SECTORS[case]
    sector = weights(5,tail,True)[index]
    if raised_spin:
        if not sector[4]:
            raise ValueError("spin is already at its maximum in this degree")
        sector = (*sector[:3],sector[3]+1,sector[4]-1)
    full = all_products(3 if group=="SU3" else 4,(*sector[3:],*sector[:3]),length,
                        include_odd_powers=group=="SU3")
    recipes = sorted({recipe for group_rows in full.values() for recipe in group_rows})
    factors = sorted({g for row in recipes for g in row})
    if group == "Sp2":
        for g in factors:
            check_factor(g)
    identifiers = {}
    keep = []
    kept_rows = []
    source_digest = hashlib.sha256()
    cursor = iter(recipes)
    while batch := list(islice(cursor,256)):
        rows = []
        for recipe in batch:
            source_digest.update(json.dumps(recipe,separators=(",",":")).encode()+b"\n")
            polynomial = coordinates(recipe,group)
            rows.append({identifiers.setdefault(k,len(identifiers)):v for k,v in polynomial.items() if v})
        all_rows = kept_rows+rows
        if not identifiers:
            continue
        dense = matrix(all_rows,len(identifiers),PRIME)
        selected = pivots(nmod_mat(dense.T.tolist(),PRIME))
        assert selected[:len(keep)] == list(range(len(keep)))
        labelled_rows = keep+batch
        keep = [labelled_rows[i] for i in selected]
        kept_rows = [all_rows[i] for i in selected]
    rank = len(keep)
    dense = matrix(kept_rows,len(identifiers),PRIME)
    cols = pivots(nmod_mat(dense.tolist(),PRIME)) if rank else []
    assert len(cols) == rank
    det = int(nmod_mat(dense[:,cols].tolist(),PRIME).det()) if rank else 1
    assert det
    labels = list(identifiers)
    checks = dict(group=group,case=case,sector=sector,index=index,length=length,prime=PRIME,
                  recipes=len(recipes),recipes_sha256=source_digest.hexdigest(),factors=len(factors),
                  monomials=len(identifiers),rank=rank,det=det,
                  selected_recipes=keep,selected_monomials=[labels[k] for k in cols])
    if raised_spin:
        checks["raised_spin"] = True
    return checks


def replay(checks):
    group = checks["group"]
    recipes = [tuple(map(tuple,r)) for r in checks["selected_recipes"]]
    keys = [tuple(tuple(map(tuple,slot)) for slot in key) for key in checks["selected_monomials"]]
    assert len(set(recipes)) == len(set(keys)) == len(keys) == len(recipes) == checks["rank"]
    positions = {k:i for i,k in enumerate(keys)}
    values = [{positions[k]:v for k,v in coordinates(row,group).items() if k in positions} for row in recipes]
    det = int(nmod_mat(matrix(values,len(keys),checks["prime"]).tolist(),checks["prime"]).det()) if keys else 1
    assert det and det == checks["det"]
    return det


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",choices=SECTORS,required=True)
    parser.add_argument("--index",type=int,choices=range(5),required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--replay",action="store_true")
    parser.add_argument("--raised-spin",action="store_true")
    args=parser.parse_args()
    if not __debug__:
        parser.error("verification requires enabled assertions")
    if args.replay:
        data=json.loads(args.out.read_bytes())
        assert data["case"]==args.case and data["index"]==args.index and data.get("raised_spin",False)==args.raised_spin
    else:
        data=record(args.case,args.index,args.raised_spin)
        raw=json.dumps(data,separators=(",",":"))+"\n"
        temporary=args.out.with_suffix(".part")
        temporary.write_text(raw)
        temporary.replace(args.out)
    print(data["case"],data["index"],data["recipes"],data["rank"],replay(data),flush=True)
