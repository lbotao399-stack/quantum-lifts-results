"""Lift a restricted SU(3) first correction through the square of the off-diagonal ghost ideal."""
import argparse
from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction
import gzip
import json
from pathlib import Path

from su3_subregular_abelian import solve_group
from su3_subregular_correction import accumulate
from su3_subregular_differential import qmonomial, variation
from su3_subregular_image import ARCHIVES, FIELDS, ROOT_START, expand, source, write


def root_degree(word):
    return len(word)-bisect_left(word, ROOT_START)


def layer_q(word):
    assert root_degree(word)==1
    return {w:n for w,n in qmonomial(word).items() if root_degree(w)==1}


def decode(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def run(level, archive, leading, recipe, build):
    data=decode(archive)
    assert data["level"]==level and data["degree"]==15 and data["linear_roots"]
    assert data["full_image_sha256"]==ARCHIVES[level][3] and data["full_image_words"]==ARCHIVES[level][1]
    rhs={bytes.fromhex(w):Fraction(-data["sign"]*c,data["denominator"]) for w,c in data["coefficients"]}
    assert all(root_degree(w)<=1 for w in rhs)
    initial=decode(leading)
    for encoded,c in initial:
        parent=bytes.fromhex(encoded)
        assert root_degree(parent)==0
        img={w:n for w,n in qmonomial(parent).items() if root_degree(w)<=1}
        accumulate(rhs,img,-Fraction(c))
    assert all(root_degree(w)==1 for w in rhs)
    original=len(rhs)
    groups=defaultdict(dict)
    inverse=defaultdict(list)
    for k in range(128):
        for w in variation(k):
            if root_degree(w)==int(k>=ROOT_START):
                assert len(w)==2 and FIELDS[k//16] in "urfg"
                inverse[w].append(k)
    for w,c in rhs.items():
        spectators=bytes(k for k in w if FIELDS[k//16] in "pqd")
        groups[spectators,len(w),w[-1]//16][w]=c
    if build:
        solution={}
        for label,group in groups.items():
            try:
                new=solve_group(group,inverse,layer_q)
            except ValueError as exc:
                raise ValueError(f"fixed leading primitive fails in component {label!r}") from exc
            if solution.keys() & new.keys():
                raise ArithmeticError("overlapping spectator components")
            solution.update(new)
        write(recipe,[[w.hex(),str(c)] for w,c in sorted(solution.items())])
    solution=decode(recipe)
    for encoded,c in solution:
        accumulate(rhs,layer_q(bytes.fromhex(encoded)),-Fraction(c))
    if rhs:
        raise ArithmeticError("nonzero residual on root-linear replay")
    return dict(level=level,full_image_words=data["full_image_words"],image_by_degree=data["nonzero_by_degree"],
                leading_terms=len(initial),linear_residual_terms=original,spectator_systems=len(groups),
                additional_terms=len(solution),remaining_terms=0)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level",choices=(32,34),required=True,type=int)
    p.add_argument("--image",required=True,type=Path)
    p.add_argument("--expand",action="store_true")
    p.add_argument("--source",type=Path)
    p.add_argument("--workers",type=int,default=2)
    p.add_argument("--leading",type=Path)
    p.add_argument("--recipe",type=Path)
    p.add_argument("--build",action="store_true")
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    if not __debug__ or args.workers<1 or not (args.expand or args.out):
        p.error("expand an image or supply a replay output, with assertions enabled")
    if args.out and (not args.leading or not args.recipe):
        p.error("a replay requires a leading primitive and a root-linear recipe")
    if args.expand:
        packet=source(args.level,args.source,15)
        packet["linear_roots"]=True
        data=expand(packet,args.workers)
        write(args.image,data)
        print("image",data["nonzero_by_degree"],flush=True)
    if args.out:
        outcome=run(args.level,args.image,args.leading,args.recipe,args.build)
        args.out.write_text(json.dumps(outcome,indent=2)+"\n")
        print(outcome,flush=True)
