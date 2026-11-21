"""Solve and replay the SU(3) scalar-chart equation in descending number of d jets."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from su3_subregular_correction import accumulate
from su3_subregular_differential import qmonomial, variation
from su3_subregular_image import ARCHIVES, write
from su3_subregular_root_coupled import d_variables, read_coefficients, solve
from su3_subregular_root_linear import decode
from su3_subregular_sparse import solve as modular_solve


@lru_cache(maxsize=60000)
def same_d(word):
    spectators = d_variables(word)
    return {w:n for w,n in qmonomial(word).items() if d_variables(w)==spectators}


def native_solution(engine, folder, rhs, depth):
    # Q0(d)=fg lowers the number of d jets. All the other variations preserve it.
    rules = [(k,term,n) for k in range(128) if k//16 != 5 for term,n in variation(k).items()]
    recipe = folder/"primitive.json.gz"
    result = modular_solve(engine, folder, rhs, recipe, None, depth, mixed_rules=rules, image_of=same_d)
    return dict(kind="primitive", coefficients=decode(recipe), systems=result["systems"])


def run(image, work, last, replay, engine=None, depth=1):
    data = decode(image)
    level = data["level"]
    assert data["degree"]==15 and data["root_limit"]==4
    assert data["sign"] == (1 if level==32 else -1)
    assert (data["full_image_words"],data["denominator"],data["full_image_sha256"])==ARCHIVES[level][1:]
    denominator = data["denominator"]
    residual = {bytes.fromhex(w):Fraction(-data["sign"]*c,denominator) for w,c in data["coefficients"]}
    original = len(residual)
    del data
    work.mkdir(parents=True,exist_ok=True)
    provenance = hashlib.sha256(image.read_bytes()).hexdigest()
    rows = []
    terms = Counter()
    for degree in range(3,last-1,-1):
        assert all(len(d_variables(w))<=degree for w in residual)
        groups = defaultdict(dict)
        for w,c in residual.items():
            s = d_variables(w)
            if len(s)==degree:
                groups[s][w]=c
        for s,rhs in sorted(groups.items(), key=lambda pair:(len(pair[1]),pair[0])):
            path=work/((s.hex() or "unit")+".json.gz")
            if path.exists() or replay:
                record = decode(path)
                assert record["full_image_file_sha256"]==provenance and record["spectator"]==s.hex()
                solution = record["solution"]
            else:
                if engine:
                    solution = native_solution(engine,work/".systems"/(s.hex() or "unit"),rhs,depth)
                else:
                    solution = solve(rhs, same_d, full=True)
                if solution["kind"] != "primitive":
                    raise ValueError(f"The higher-d choices require a joint solve at {s.hex()}")
                record=dict(full_image_file_sha256=provenance, spectator=s.hex(), solution=solution)
            recipe = read_coefficients(solution["coefficients"])
            assert all(d_variables(w)==s for w in recipe)
            check = dict(rhs)
            for w,c in recipe.items():
                full = qmonomial(w)
                accumulate(check, {m:n for m,n in full.items() if d_variables(m)==s}, -c)
                accumulate(residual, full, -c)
            if check:
                raise ArithmeticError("Nonzero rational residual on the selected d coefficient")
            if not path.exists():
                temporary = path.with_suffix(".part.gz")
                write(temporary,record)
                temporary.replace(path)
            terms[degree]+=len(recipe)
            rows.append(dict(spectator=s.hex(),input=len(rhs),primitive=len(recipe)))
            print("level",level,"d",degree,s.hex(),"input",len(rhs),"primitive",len(recipe),flush=True)
        assert all(len(d_variables(w))<degree for w in residual)
        print("through d",degree,"remaining",len(residual),flush=True)
    return dict(level=level,full_image_file_sha256=provenance,original_terms=original,
                through_d_count=last,primitive_by_d_count=dict(terms),remaining_terms=len(residual),components=rows)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--image",type=Path,required=True)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--through",type=int,choices=range(4),default=0)
    p.add_argument("--replay",action="store_true")
    p.add_argument("--engine",type=Path,help="Compiled subregular_sparse.cpp for mixed Taylor degrees")
    p.add_argument("--depth",type=int,default=1,help="Incidence rounds for the native solve")
    args=p.parse_args()
    if not __debug__ or args.depth<1 or (args.replay and args.engine):
        p.error("Enable assertions and select either a solve or a replay")
    result=run(args.image,args.work,args.through,args.replay,args.engine,args.depth)
    args.out.write_text(json.dumps(result,indent=2)+"\n")
