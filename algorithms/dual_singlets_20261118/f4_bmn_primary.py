"""Classical bounds in the F4 level-24 derivative-free source selected by six commutant charges."""
import argparse
import math
from collections import Counter
from pathlib import Path

import numpy as np
from flint import nmod_mat

import mixed_incoming as m
import scalar_neighbors as scalar
from f4_colour import frame
from f4_dressed_jordan import check_short_arms, source, variation
from f4_formal_cycles import coefficient_matrix,probe
from f4_joined_boundaries import evaluate,native
from f4_singlet_characters import singlets
from fortuity.full.letters import Q_letter
from mixed33_exact import read

LETTERS=tuple((0,0,t) for t in ((0,),(1,),(2,),(1,2),(0,2),(0,1),(0,1,2)))
SIGNS=(1,1,1,1,-1,1,-1)
WEIGHT=(0,0,4,4,4)


def bmn_convention():
    for i,l in enumerate(LETTERS):
        left=Counter({(a,b):SIGNS[i]*n for n,a,b in Q_letter(l)})
        right=Counter({(LETTERS[a],LETTERS[b]):n*SIGNS[a]*SIGNS[b] for n,a,b in scalar.Q.get(i,())})
        assert left==right


def translated(expr):
    out=Counter()
    for w,c in expr.items():
        coeff=c*math.prod(SIGNS[i] for tr in w for i in tr)
        sg,v=m.so_canon(tuple(tuple(LETTERS[i] for i in tr) for tr in w))
        if sg:
            out[v]+=sg*coeff
    return {w:c for w,c in out.items() if c}


def run(work,witness,tasks,workers):
    work.mkdir(exist_ok=True,parents=True)
    bmn_convention()
    checked=check_short_arms()
    data=read(witness)
    p=data["prime"]
    results=[]
    trace_lib=m.library(work,27,p)
    jordan=native(work,p,"jordan")
    for length in (11,10):
        dim=singlets("F4",WEIGHT,length+1)
        assert len(dim)==1
        degree=next(iter(dim))
        definitions=list(dict.fromkeys(scalar.decode(rec)
            for group in ("SO_odd","Sp") for rec in data["maps"][f"{group}:{length}"]["source_contractions"]))
        exprs=[translated(scalar.expression(r)) for r in definitions]
        assert all(not m.q(e) for e in exprs)
        words,rows=m.flatten(exprs)
        stored,frames=probe(rows,words,dim,p,trace_lib,workers,715180+length*73)
        base=m.minor(stored,p)
        counts=[base["rank"]]
        alphabet={l for w in words for tr in w for l in tr}
        ff=[frame(alphabet,x["points"],x["k"],degree,p,x["seed"]) for x in frames]
        mats=[stored]
        seeds=[]
        start=871161+1037*length
        total_graphs=[]
        # A squared Gram determinant of three independent Cartan scalars survives
        # at length 12. Results without this bound still use independent minors.
        for i in range(start,start+tasks,64):
            recipes=[]
            for seed in range(i,min(i+64,start+tasks)):
                row=source(0,4,12-length,seed)
                if not row:
                    continue
                assert all(not variation(row,k) for k in (1,2,3))
                image=variation(row,0)
                if image:
                    assert not variation(image,0)
                    recipes.append(image);seeds.append(seed)
            if not recipes:
                continue
            graphs=sorted({g for e in recipes for g in e});total_graphs.append(len(graphs))
            assert all(m.nferm(g)==degree and all(l in alphabet for tr in g for l in tr) for g in graphs)
            ids={g:i for i,g in enumerate(graphs)}
            co=coefficient_matrix([[(ids[g],n) for g,n in e.items()] for e in recipes],len(graphs),p)
            mats.append(np.column_stack([np.asarray(co@evaluate(graphs,f,work,workers,jordan,"jordan"))%p for f in ff]))
            minor=m.minor(np.vstack(mats),p)
            counts.append(minor["rank"])
            print(length,"Jordan sources",len(seeds),"rank",minor["rank"],"ordinary",base["rank"],flush=True)
            if length==11 and minor["rank"]==dim[degree]-1:
                break
            if length==10 and minor["rank"]+results[0]["minor"]["rank"]==dim[degree]:
                break
        if len(mats)==1:
            minor=base
        results.append(dict(length=length,rank_history=counts,ordinary=base,minor=minor,
                            ordinary_recipes=definitions,jordan_seeds=seeds,frames=frames,graphs=total_graphs))
        dims=[sum(singlets("F4",WEIGHT,l).values()) for l in (10,11,12)]
        upper=dims[1]-sum(x["minor"]["rank"] for x in results) if len(results)==2 else None
        assert upper is None or upper>=0
        m.write(work/"primary.json",dict(prime=p,weight=WEIGHT,dimensions=dims,short_arm_checks=checked,
                 letter_action_checks=len(LETTERS),classical_cohomology_upper=upper,records=results))
    return upper


if __name__=="__main__":
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--work",type=Path,required=True)
    a.add_argument("--witness",type=Path,default=Path(__file__).with_name("scalar_minors.json.gz"))
    a.add_argument("--tasks",type=int,default=640)
    a.add_argument("--workers",type=int,default=16)
    args=a.parse_args();assert __debug__
    print("classical bound",run(**vars(args)),flush=True)
