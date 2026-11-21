"""Full first equation on the subregular chart with both off-diagonal ghosts f,g set to zero."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import partial
import gzip
import json
from pathlib import Path

from su3_subregular_differential import FIELDS, ODD, qmonomial, variation
from su3_subregular_image import source, expand, write
from su3_subregular_correction import accumulate, quadratic_parents, reduce_row


def survives(word, vanish):
    return all(FIELDS[k//16] not in vanish for k in word)


def projected_q(word, vanish):
    return {w:n for w,n in qmonomial(word).items() if survives(w, vanish)}


def solve_group(rhs, inverse, differential):
    # Every term in a connected system has the same jets of the closed spectators.
    reached, parents = set(rhs), set()
    frontier=set(rhs)
    while frontier:
        incident={p for w in frontier for p in quadratic_parents(w,inverse)}-parents
        parents.update(incident)
        incoming={w for p in incident for w in differential(p)}
        frontier=incoming-reached
        reached.update(frontier)
    pivots={}
    recipes=[]
    for p in sorted(parents):
        image=differential(p)
        steps=reduce_row(image,pivots)
        if image:
            first=min(image)
            factor=image[first]
            pivots[first]=len(recipes), {w:Fraction(n,factor) for w,n in image.items()}
            recipes.append((p,factor,steps))
    reduced=dict(rhs)
    steps=reduce_row(reduced,pivots)
    if reduced:
        raise ValueError(f"nonzero quotient in a complete spectator component: {len(reduced)} terms")
    coefficients=Counter(dict(steps))
    primitive={}
    for i in reversed(range(len(recipes))):
        if coefficients[i]:
            p,den,subs=recipes[i]
            c=Fraction(coefficients[i],den)
            primitive[p]=c
            for j,n in subs:
                coefficients[j]-=n*c
    return primitive


def run(level, image, primitive, build, workers, vanish="dfg", input_archive=None):
    # There are seven B factors of degree at most one and four A factors.
    # With d=0, A has degree at most one. With d present, A=dB+a diag(1,1,-2).
    degree = 11 if vanish == "dfg" else 15
    if build:
        packet=source(level,input_archive,degree,vanish=vanish)
        result=expand(packet,workers)
        write(image,result)
    result=json.loads(gzip.decompress(image.read_bytes()))
    assert result["level"]==level and result["degree"]==degree and result["vanishing_fields"]==vanish
    rhs={bytes.fromhex(w):Fraction(-result["sign"]*n,result["denominator"]) for w,n in result["coefficients"]}
    assert all(survives(w, vanish) for w in rhs)
    inverse=defaultdict(list)
    for k in range(128):
        terms={w:n for w,n in variation(k).items() if survives(w, vanish)}
        if not survives(bytes([k]), vanish):
            assert not terms
        else:
            for w in terms:
                assert len(w)==2 and FIELDS[k//16] in "ur"
                inverse[w].append(k)
    groups=defaultdict(dict)
    for w,n in rhs.items():
        frozen=bytes(k for k in w if FIELDS[k//16] in "pqd")
        groups[frozen, len(w)][w]=n
    if build:
        solution={}
        differential=partial(projected_q, vanish=vanish)
        for group in groups.values():
            solution.update(solve_group(group,inverse,differential))
        write(primitive,[[w.hex(),str(c)] for w,c in sorted(solution.items())])
    solution=json.loads(gzip.decompress(primitive.read_bytes()))
    for word,c in solution:
        w=bytes.fromhex(word)
        assert survives(w, vanish)
        accumulate(rhs,projected_q(w, vanish),-Fraction(c))
    assert not rhs
    return dict(level=level,full_image_words=result["full_image_words"],monomials=result["nonzero_by_degree"],
                spectator_groups=len(groups),primitive_terms=len(solution),remaining_terms=0)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level",type=int,choices=(32,34),required=True)
    p.add_argument("--image",type=Path,required=True)
    p.add_argument("--primitive",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--build",action="store_true")
    p.add_argument("--workers",type=int,default=1)
    p.add_argument("--vanish",choices=("dfg","fg"),default="dfg")
    p.add_argument("--source",type=Path)
    args=p.parse_args()
    if not __debug__ or args.workers < 1:
        p.error("verification requires assertions and a positive worker count")
    result=run(args.level,args.image,args.primitive,args.build,args.workers,args.vanish,args.source)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
