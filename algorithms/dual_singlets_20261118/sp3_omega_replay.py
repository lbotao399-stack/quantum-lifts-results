"""Verify the rational omega cycle on new integer jets and project its lowered weight."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import mixed_incoming as m
from dual_charge_bands import coordinates
from full_extremal import project
from so8_level24_descendant import raising_images, substitute
from sp3_omega_actions import centre_basis, load_page, signed_frame
from sp3_omega_substitutions import STEPS
from trace_descent import evaluate


def pair(polynomial, data, library, workers):
    weight = tuple(data["weight"])
    letters = sorted({l for content in m.contents(weight,data["length"]) for l in content})
    values = np.concatenate([coordinates([polynomial],int(degree),letters,
                frame["samples"],6,"Sp",library,workers)[0] for degree,frame in data["centre_frames"].items()])
    functional = np.asarray(data["representative"]["target_functional"],dtype=np.int64)
    assert len(values)==len(functional) and len(values)*(m.P-1)**2<2**63
    return int(values@functional % m.P)


def run(page_path, cycle_path, lower_path, out, workers, seed):
    page = load_page(page_path)
    cycle, lower = [json.loads(gzip.decompress(p.read_bytes())) for p in (cycle_path,lower_path)]
    assert cycle["page_sha256"]==hashlib.sha256(page_path.read_bytes()).hexdigest()
    assert (lower["weight"],lower["length"],lower["group"],lower["colours"])==([1,1,5,2,2],9,"Sp",6)
    rows = centre_basis(page)
    coefficients = [(i,Fraction(q)) for i,q in cycle["cycle_basis_coefficients"]]
    common = math.lcm(*(c.denominator for _,c in coefficients))
    numerator=Counter()
    for i,c in coefficients:
        scale=int(c*common)
        for w,n in rows[i].items():
            numerator[w]+=scale*n
    numerator={w:c for w,c in numerator.items() if c}
    assert not any(any(r.values()) for r in raising_images(numerator))
    out.mkdir(parents=True,exist_ok=True)
    library=m.library(out,6,m.P)
    frame=signed_frame(tuple(page["weight"]),10,seed)
    proposal=m.q(numerator)
    closed=evaluate([proposal],frame,6,"Sp")[0]
    assert not any(closed)
    assert pair(numerator,page,library,workers) == common % m.P
    original=numerator
    for spin,flavour in STEPS:
        numerator=substitute(numerator,spin,flavour,False)
    raw_raising=[sum(bool(v) for v in r.values()) for r in raising_images(numerator)]
    quotient=project(numerator,m.sector(tuple(lower["weight"]),9))
    assert not any(any(r.values()) for r in raising_images(quotient))
    pvalue=pair(quotient,lower,library,workers)*pow(common,-1,m.P)%m.P
    print("lowered",len(numerator),"raising supports",raw_raising,"projected",len(quotient),
          "pairing",pvalue,flush=True)
    for spin,flavour in reversed(STEPS):
        quotient=substitute(quotient,spin,flavour,True)
    result=project(quotient,m.sector(tuple(page["weight"]),9))
    returned=pair(result,page,library,workers)*pow(common,-1,m.P)%m.P
    first=min(original)
    coefficient=Fraction(result.get(first,0),original[first])
    same=all(result.get(w,0)==coefficient*original.get(w,0) for w in result.keys()|original.keys())
    record=dict(page_sha256=cycle["page_sha256"],cycle_sha256=hashlib.sha256(cycle_path.read_bytes()).hexdigest(),
                lower_sha256=hashlib.sha256(lower_path.read_bytes()).hexdigest(),common_denominator=common,
                cycle_terms=len(original),formal_variation_terms=len(proposal),signed_seed=seed,
                signed_frames=frame,verified_closure_coordinates=len(closed),closure_nonzero=0,
                original_pairing_mod_prime=1,raw_lowered_terms=len(numerator),raw_raising_terms=raw_raising,
                projected_lower_pairing_mod_prime=pvalue,projected_return_pairing_mod_prime=returned,
                exact_formal_return_coefficient=str(coefficient) if same else None,prime=m.P)
    m.write(out/"cycle_replay.json.gz",record)
    print("closed coordinates",len(closed),"projected return",record["exact_formal_return_coefficient"],
          "residue",returned,flush=True)
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for opt in ("page","cycle","lower","out"):
        parser.add_argument("--"+opt,type=Path,required=True)
    parser.add_argument("--workers",type=int,default=4)
    parser.add_argument("--seed",type=int,default=600911)
    args=parser.parse_args()
    assert __debug__
    run(args.page,args.cycle,args.lower,args.out,args.workers,args.seed)
