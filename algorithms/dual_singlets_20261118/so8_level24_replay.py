"""Replay the rational SO(8) corrections from trace words on new integer colour frames."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import mixed_incoming as m
from so8_level24_generator import even_class
from trace_descent import action, basis, dimensions, evaluate, highest_weights, settings


def selected_actions(grade, band, chain, operator):
    if grade == 0:
        row = even_class(operator)
        return m.q(row), action(row,1,8,"SO")
    record = band["frames"][grade]
    formal = highest_weights(tuple(record["weight"]),record["length"]-1)
    parents = [formal[int(degree)][i] for degree,frame in record["source"].items() for i in frame["rows"]]
    assert len(parents)==band["source_dimensions"][grade]
    classical, quantum = Counter(), Counter()
    for index, coefficient in chain["corrections"][grade-1]["classical_rows"]:
        row=parents[index]
        first=action(row,1,8,"SO")
        scale=math.lcm(*(Fraction(v).denominator for v in first.values()))
        n=Fraction(coefficient)*scale
        for w,c in m.q(row).items(): classical[w]+=n*c
        for w,c in first.items(): quantum[w]+=n*c
    return tuple({w:c for w,c in image.items() if c} for image in (classical,quantum))


def worker(task):
    grade,band,chain,op,seed=task
    settings()
    _,advance = selected_actions(grade-1,band,chain,op) if grade else ({},{})
    diagonal,_ = selected_actions(grade,band,chain,op)
    residual=Counter(diagonal)
    residual.update(advance)
    row={w:n for w,n in residual.items() if n}
    denominator=math.lcm(*(Fraction(n).denominator for n in row.values()))
    integer={w:int(n*denominator) for w,n in row.items()}
    record=band["frames"][grade]
    weight,length=tuple(record["weight"]),record["length"]
    target=highest_weights(weight,length)
    dims=dimensions(weight,length,8,"SO")
    alphabet=sorted({l for content in m.contents(weight,length) for l in content})
    _,frames=basis(target,dims,"SO",8,seed+1009*grade,alphabet)
    values=evaluate([integer],frames,8,"SO")[0]
    assert len(values)==sum(dims.values()) and not any(values)
    result=dict(grade=grade,universal_residual_terms=len(row),common_denominator=denominator,
                exact_coordinates=len(values),nonzero_coordinates=sum(bool(v) for v in values),
                character_dimensions=dims,frame_determinants={k:v["determinant_mod_prime"] for k,v in frames.items()})
    print(result,flush=True)
    return result


def run(operator,band_path,chain_path,workers,seed):
    band=json.loads(band_path.read_bytes())
    chain=json.loads(chain_path.read_bytes())
    assert hashlib.sha256(band_path.read_bytes()).hexdigest()==chain["band_sha256"]
    assert hashlib.sha256(operator.read_bytes()).hexdigest()==chain["operator_sha256"]
    assert [s["grade"] for s in chain["corrections"]]==[1,2,3,4]
    tasks=[(grade,band,chain,operator,seed) for grade in range(5)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        grades=list(pool.map(worker,tasks))
    return dict(seed=seed,chain_sha256=hashlib.sha256(chain_path.read_bytes()).hexdigest(),
                frame_prime=m.P,grades=grades)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    folder=Path(__file__).with_name("so8_level24")
    parser.add_argument("--operator",type=Path,default=folder/"even_bmn.json.gz")
    parser.add_argument("--band",type=Path,default=folder/"band1_4_4_4_s0_so8.json")
    parser.add_argument("--chain",type=Path,default=folder/"chain.json")
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--seed",type=int,default=570331)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert __debug__
    m.write(args.out,run(args.operator,args.band,args.chain,args.workers,args.seed))
