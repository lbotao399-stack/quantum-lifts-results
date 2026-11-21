"""Bound odd scalar cohomology by restricting full-group adjoint trace boundaries to F4."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from bmn_primary_index import cochains
import mixed_incoming as m
from f4_colour import frame
from f4_adjoint_traces import adjoint_frame


def run(weights, work, prime, seed, points):
    if any(n%2==0 for n in map(sum,weights)):
        raise ValueError("Odd total scalar degree is required")
    letters = [(0,0,(i,)) for i in range(3)]
    work.mkdir(parents=True,exist_ok=True)
    frames = {27:frame(letters,points,1,0,prime,seed),52:adjoint_frame(letters,points,1,0,prime,seed)}
    libraries = {k:m.library(work,k,prime) for k in frames}
    maximum = max(map(sum,weights))
    cochain = {g:cochains(maximum,g)[1] for g in ("F4","E6","E7","E8")}
    rows = []
    for weight in weights:
        length = sum(weight)
        blocks = m.words_by_shape(length,(0,0,*weight))
        words, highest = m.flatten([e for _,block in sorted(blocks.items()) for e in m.singlets(block)])
        if not highest:
            continue
        traces = sorted({t for w in words for t in w})
        index = {t:i for i,t in enumerate(traces)}
        primitive = {k:m.values([(t,) for t in traces],v,libraries[k]) for k,v in frames.items()}
        result = {}
        for group,copies in (("F4",0),("E6",1),("E7",3),("E8",7)):
            # F4 uses its 26; a trivial line in the 27 annihilates nonempty words.
            tr = primitive[27] if group == "F4" else (primitive[52]+copies*primitive[27])%prime
            coeff = np.empty((len(words),points),dtype=np.int64)
            for j,w in enumerate(words):
                value=np.ones(points,dtype=np.int64)
                for t in w:
                    value=value*tr[index[t]]%prime
                coeff[j]=value
            evaluated = np.zeros((len(highest),points),dtype=np.int64)
            for j,row in enumerate(highest):
                for k,n in row:
                    evaluated[j]=(evaluated[j]+(n%prime)*coeff[k])%prime
            witness=m.minor(evaluated,prime)
            dimensions=cochain[group][length].get(weight,[0,0] if group=="E6" else [0])
            assert witness["rank"]<=dimensions[0]
            result[group]={"dimension":dimensions[0], "boundary_minor":witness,
                           "minor_matrix":evaluated[np.ix_(witness["rows"],witness["columns"])].tolist(),
                           "cohomology_upper":dimensions[0]-witness["rank"]}
        rows.append({"weight":weight,"length":length,"groups":result})
        print(weight, {g:(r['dimension'],r['boundary_minor']['rank']) for g,r in result.items()},flush=True)
    return dict(prime=prime,seed=seed,points=points,color_component_E6="outer even",results=rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--prime",type=int,choices=(1000003,1000033),default=1000003)
    p.add_argument("--seed",type=int,default=334057)
    p.add_argument("--points",type=int,default=22)
    p.add_argument("--weight",action="append",nargs=3,type=int,required=True)
    a=p.parse_args()
    data=run([tuple(w) for w in a.weight],a.work,a.prime,a.seed,a.points)
    a.out.write_text(json.dumps(data,separators=(",",":"))+"\n")
